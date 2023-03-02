# -*- coding: utf-8 -*-
"""Core implementation of the W-Phase inversion."""
import logging
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from functools import partial
from time import perf_counter
from typing import Any, List, Mapping, NamedTuple, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from numpy.linalg import lstsq
from numpy.typing import ArrayLike
from obspy import Trace
from obspy.core import Stream, UTCDateTime
from scipy.interpolate import interp1d
from scipy.ndimage import convolve1d
from scipy.optimize import OptimizeResult, minimize

from wphase.wputils import SimpleTimer

try:
    from scipy.signal.windows import triang
except ImportError:
    from scipy.signal import triang

try:
    from obspy.geodetics import locations2degrees
except ImportError:
    from obspy.core.util.geodetics import locations2degrees

try:
    from obspy.geodetics.base import gps2dist_azimuth
except ImportError:
    from obspy.core.util.geodetics import gps2DistAzimuth as gps2dist_azimuth

try:
    from obspy.signal.invsim import paz_2_amplitude_value_of_freq_resp
except ImportError:
    from obspy.signal.invsim import paz2AmpValueOfFreqResp as paz_2_amplitude_value_of_freq_resp

from wphase import settings

from .bandpass import bandpassfilter
from .exceptions import (DeconvolutionError, InversionError,
                         UnknownTransferFunction)
from .greens import AbsoluteComponent, GreensFunctions
from .model import (ChannelMetadata, Event, OL1Result, OL2Result, OL3Result,
                    PreliminaryMagnitudeFit, TimeDelayMisfits, WPhaseResult)
from .seismoutils import azimuthal_gap, ltrim, pad1d, rot_12_NE

logger = logging.getLogger(__name__)

@dataclass
class ChannelData:
    """We just use this to provide type hints for rows of our `data` frame."""
    name: str
    metadata: ChannelMetadata
    location: str
    component: AbsoluteComponent
    raw_trace: Trace
    trace: Trace
    distance: float
    azimuth: float
    ptime: float
    start_time: UTCDateTime
    end_time: UTCDateTime


def build_dataframe(
    stream: Stream,
    event: Event,
    metadata: Mapping[str, ChannelMetadata],
    ptimes: Mapping[str, float]
) -> pd.DataFrame:
    """Combine waveform data and our extra associated data into a single dataframe."""
    #In case of multiple locations we favor
    # '00'. This may be improved if we select the one with the longer period
    # sensor.
    st_sel = stream.select(location = '00')
    st_sel += stream.select(location = '--')
    st_sel += stream.select(location = '') #Check this.

    # We don't really care exactly what ordering obspy chooses
    # here; we just want a deterministic ordering so we can compare runs.
    st_sel.sort()

    logger.info(
        'Initial number of traces: using %d with default locations (of %d total)',
        len(st_sel), len(stream)
    )

    comps = settings.WAVEFORM_COMPONENTS
    if comps:
        st_sel = Stream([tr for tr in st_sel if tr.stats.component in comps])
        logger.info('%d traces remaining after restricting to %s', len(st_sel), comps)

    # rotate the horizontal components to geographical north or east
    st_sel = rot_12_NE(st_sel, metadata)

    # Convert obspy Stream to a dataframe so we can attach our extra data to each channel
    def trace_to_record(trace: Trace):
        return {
            "id": trace.id,
            "location": trace.stats.location,
            "component": trace.stats.component,
            "raw_trace": trace,
            "metadata": metadata.get(trace.id),
        }
    data = pd.DataFrame.from_records(map(trace_to_record, st_sel), index="id")

    def epicentral_distance(row):
        return locations2degrees(event.latitude, event.longitude, row.metadata.latitude, row.metadata.longitude)
    data["distance"] = data.apply(epicentral_distance, axis=1)

    def azimuth_from_north(row):
        return gps2dist_azimuth(event.latitude, event.longitude, row.metadata.latitude, row.metadata.longitude)[1]
    data["azimuth"] = data.apply(azimuth_from_north, axis=1)

    gap = azimuthal_gap(data.azimuth)
    if gap > settings.MAXIMUM_AZIMUTHAL_GAP:
        raise InversionError("Lack of azimuthal coverage (%.0f° > %.0f°). Aborting."
                             % (gap, settings.MAXIMUM_AZIMUTHAL_GAP))

    def travel_time(row):
        return ptimes.get(row.name)
    data["ptime"] = data.apply(travel_time, axis=1)
    data = data[data["ptime"].notnull()].copy()

    origin_time = UTCDateTime(event.time)
    data["start_time"] = data["ptime"].map(lambda t: origin_time + t)
    data["end_time"] = data["start_time"] + data["distance"] * settings.WPHASE_CUTOFF

    return data

_Vpaz2freq = np.vectorize(paz_2_amplitude_value_of_freq_resp)
def amplitude_at_frequency(metadata: ChannelMetadata, frequency: np.ndarray):
    """Given metadata for a seismometer channel, return the amplitude response
    at a given frequency."""
    if metadata.transfer_function == "B":
        frequency = frequency / (2*np.pi)
    elif metadata.transfer_function != "A":
        raise UnknownTransferFunction(metadata.transfer_function)
    return _Vpaz2freq({
        "poles": metadata.poles,
        "zeros": metadata.zeros,
        "gain": metadata.gain,
    }, frequency)

def preprocess_channel_waveforms(row: ChannelData, frequencies: np.ndarray, Ta: float, Tb: float) -> Optional[Trace]:
    """Deconvolve + filter the waveform data for a single channel."""
    trid = row.name

    try:
        response_amplitude = amplitude_at_frequency(row.metadata, frequencies)
    except UnknownTransferFunction:
        logger.warning("Unknown transfer function. Skipping %s", trid)
        return None

    # Fitting the instrument response and getting coefficients
    response_coefficients = fit_instrument_response(
        row.metadata.sensitivity, frequencies, response_amplitude
    )
    if response_coefficients is None:
        logger.warning("Couldn't fit instrument response coefficients. Skipping %s", trid)
        return None
    else:
        om0, h, G = response_coefficients

    omega = frequencies*2.*np.pi
    amplitude_from_coeff = np.abs(omega*omega / (omega*omega + 2j*h*om0*omega - om0*om0))

    # Misfit is the relative L2 norm of the difference, expressed as a percentage
    misfit = (
        100*np.linalg.norm(response_amplitude-amplitude_from_coeff)
        / np.linalg.norm(response_amplitude)
    )

    if misfit > settings.RESPONSE_MISFIT_TOL:
        logger.warning('Bad fit for instrument response (misfit=%e). Skipping %s', misfit, trid)
        return None

    dt = row.raw_trace.stats.delta
    try:
        ret = Trace(
            deconvolve_and_filter(
                row.raw_trace,
                om0,
                h,
                G,
                dt,
                corners=4,
                baselinelen=60./dt,
                taperlen=10.,
                fmin=1./Tb,
                fmax=1./Ta
            ),
            row.raw_trace.stats,
        )

    except DeconvolutionError as e:
        logger.warning("Error deconvolving trace %s: %s", trid, str(e))
        return None

    # trim to the Wphase time window
    ret.trim(row.start_time, row.end_time)
    if len(ret) == 0:
        logger.warning("Empty trace. Skipping %s", trid)
        return None

    return ret

def waveform_preprocessor(Ta: float, Tb: float):
    T = np.linspace(Ta,Tb,500)
    freq = 1./T
    return partial(preprocess_channel_waveforms, frequencies=freq, Ta=Ta, Tb=Tb)


def reject_p2p_outliers(data: pd.DataFrame) -> pd.DataFrame:
    median_p2p = data.p2p.median(skipna=True)
    mrcoeff = settings.MEDIAN_REJECTION_COEFF
    p2pmax, p2pmin = mrcoeff[1]*median_p2p, mrcoeff[0]*median_p2p

    for trid, row in data.iterrows():
        if np.isnan(row.p2p):
            logger.warning("P2P Amp for %s is NaN! Excluding.", trid)
        elif row.p2p > p2pmax:
            logger.info("P2P Amp for %s is too big (%.2e > %.2e). Excluding.", trid, row.p2p, p2pmax)
        elif row.p2p < p2pmin:
            logger.info("P2P Amp for %s is too small (%.2e < %.2e). Excluding.", trid, row.p2p, p2pmin)

    return data[data.p2p.between(p2pmin, p2pmax)] # type: ignore


def computeOL1(data: pd.DataFrame, result: WPhaseResult):
    """Preliminary magnitude estimate based only on P2P amplitudes in W-Phase
    spectral band."""
    data = data.copy()
    data["trace"] = data.apply(waveform_preprocessor(200, 1000), axis=1)

    # Remove channels that failed waveform preprocessing:
    data = data[data.trace.notnull()].copy() # type: ignore
    data["p2p"] = data.trace.map(lambda d: d[:].max() - d[:].min())

    data = reject_p2p_outliers(data[data.component == "Z"]).sort_values("distance") # type: ignore

    if len(data) < settings.MINIMUM_STATIONS:
        raise InversionError("Lack of stations (%d < %d). Aborting."
                             % (len(data), settings.MINIMUM_STATIONS))

    logger.info("Traces accepted for preliminary magnitude calculation: {}"
                .format(len(data)))

    pre_results = preliminary_magnitude(data.p2p, data.distance, data.azimuth, data.index)

    logger.info("Preliminary t_h = %.2f" % pre_results.timescale)

    logger.info("OL1:")
    logger.info("Average amplitude:  %.1em", pre_results.average_amplitude)
    logger.info("Magnitude:          %.2f", pre_results.magnitude)
    logger.info("Strike:             %.1f°", pre_results.strike)
    logger.info("Eccentricity (b/a): %.2f", pre_results.eccentricity)

    result.OL1 = OL1Result(
        preliminary_calc_details=pre_results,
        used_traces=list(data.index),
        nstations=len(data),
        magnitude=round(pre_results.magnitude, 1),
    )


def computeNew(
    data: pd.DataFrame,
    event: Event,
    greens: GreensFunctions,
    result: WPhaseResult,
    t_h: float,
):
    ol1 = result.OL1
    if ol1 is None:
        raise ValueError("OL1 result not found, can't proceed with optimization")

    # Now that we have a preliminary magnitude, we use it to choose a more precise frequency band:
    Ta, Tb = get_corner_freqs_from_mag(ol1.preliminary_calc_details.magnitude)
    logger.info("OLNew Filter corner periods: %.1f, %1f" % (Ta, Tb))

    data = data.copy()
    data["trace"] = data.apply(waveform_preprocessor(Ta, Tb), axis=1)
    data = data[data.trace.notnull()].copy() # type: ignore
    data["p2p"] = data.trace.map(lambda d: d[:].max() - d[:].min())
    data = reject_p2p_outliers(data).sort_values("distance")

    logger.info('number of traces for OLNew calculation: %d', len(data))
    params = Point(event.latitude, event.longitude, event.depth)
    for tol in settings.MISFIT_TOL_SEQUENCE:
        _, (M, misfit, GFmatrix) = omni_minimize(
            data, greens, (Ta, Tb), t_h, t_h, params, tol=settings.OPTIMIZATION_TOLERANCE
        )
        logger.info("With %d channels, overall misfit is %.0f", len(data), misfit)
        syn = GFmatrix @ full_to_deviatoric(M)
        def misfit_of_channel(row):
            channel_syn = syn[row.start:row.end]
            return 100 * np.linalg.norm(channel_syn - row.trace) / np.linalg.norm(channel_syn)
        misfits = data.join(trace_indices(data.trace)).apply(misfit_of_channel, axis=1)
        data = data[misfits <= tol]

        logger.info("After culling to misfit <= %.0f%%: %d channels remaining", tol, len(data))
        if len(data) < settings.MINIMUM_FITTING_CHANNELS:
            msg = "Only {} channels with possibly acceptable fits. Aborting.".format(len(data))
            logger.error(msg)
            raise InversionError(msg, result=result)

    params, (M, misfit, GFmatrix) = omni_minimize(data, greens, (Ta, Tb), t_h, t_h, params)
    syn = GFmatrix @ full_to_deviatoric(M)
    trace_lengths = list(zip(data.index, data.trace.map(len)))
    result.new = make_result(
        OL3Result,
        M,
        misfit=misfit,
        depth=params.depth,
        time_delay=t_h,
        centroid=(params.latitude, params.longitude, params.depth),
        periods=(Ta, Tb),
        used_traces=list(data.index),
        moment_tensor=M,
        observed_displacements=np.concatenate(data.trace),
        synthetic_displacements=syn,
        trace_lengths=OrderedDict(trace_lengths),
        grid_search_candidates=[],
        grid_search_results=[],
    )

def computeOL2(
    data: pd.DataFrame,
    event: Event,
    greens: GreensFunctions,
    result: WPhaseResult,
):
    """Moment tensor from PDE inversion, using fixed centroid (from preliminary
    hypocenter) and searching for optimal time delay."""
    ol1 = result.OL1
    if ol1 is None:
        raise ValueError("OL1 result not found, can't proceed with OL2")

    # Now that we have a preliminary magnitude, we use it to choose a more precise frequency band:
    Ta, Tb = get_corner_freqs_from_mag(ol1.preliminary_calc_details.magnitude)
    logger.info("OL2 Filter corner periods: %.1f, %1f" % (Ta, Tb))

    data = data.copy()
    data["trace"] = data.apply(waveform_preprocessor(Ta, Tb), axis=1)

    # Remove channels that failed waveform preprocessing:
    data = data[data.trace.notnull()].copy() # type: ignore
    data["p2p"] = data.trace.map(lambda d: d[:].max() - d[:].min())
    data = reject_p2p_outliers(data).sort_values("distance")

    logger.info('number of traces for OL2: %d', len(data))
    trlen = np.asarray(data.trace.map(len))
    mrf = moment_rate_function(ol1.preliminary_calc_details.timescale, greens.delta)
    max_t_d = int(settings.MAXIMUM_TIME_DELAY * greens.delta)
    time_delays = np.arange(1., max_t_d)

    centroid = Point(
        latitude=event.latitude,
        longitude=event.longitude,
        depth=event.depth,
    )

    # We start by minimizing misfit over time delay. Since changing time delay
    # just translates the waveforms, we only need to fetch the greens
    # functions once:
    GFmatrix = get_greens_for_dataset(
        greens=greens,
        data=data,
        centroid=centroid,
        periods=(Ta, Tb),
        mrf=mrf,
        padding=max_t_d,
        t_d=max_t_d,
    )

    # This was taking *longer* to run when we used a parallel pool - the
    # overhead doesn't seem to be worth it.
    observed_displacements = np.concatenate(data.trace)
    misfits = [
        get_timedelay_misfit(t, GFmatrix, trlen, observed_displacements, max_t_d)
        for t in time_delays
    ]

    # Set t_d (time delay) and and t_h (half duration) to optimal values:
    mis_min = int(np.array(misfits).argmin())
    result.misfits = TimeDelayMisfits(array=misfits, min=mis_min)

    t_d = t_h = time_delays[mis_min]
    mrf = moment_rate_function(t_h, greens.delta)
    logger.info("Source time function, time delay: %d, %f", len(mrf), t_d)
    logger.info("revised t_h = %.2f" % t_h)

    inversion_params = dict(
        t_d=t_d,
        greens=greens,
        centroid=centroid,
        periods=(Ta, Tb),
        mrf=mrf,
    )

    #### Removing individually bad-fitting channels.
    # This iteratively removes channels with misfits outside of the acceptable
    # range defined by the setting MISFIT_TOL_SEQUENCE.
    for tol in settings.MISFIT_TOL_SEQUENCE:
        M, misfit, GFmatrix = core_inversion(data=data, **inversion_params)
        logger.info("With %d channels, overall misfit is %.0f", len(data), misfit)
        syn = GFmatrix @ full_to_deviatoric(M)
        def misfit_of_channel(row):
            channel_syn = syn[row.start:row.end]
            return 100 * np.linalg.norm(channel_syn - row.trace) / np.linalg.norm(channel_syn)
        misfits = data.join(trace_indices(data.trace)).apply(misfit_of_channel, axis=1)
        data = data[misfits <= tol]

        logger.info("After culling to misfit <= %.0f%%: %d channels remaining", tol, len(data))
        if len(data) < settings.MINIMUM_FITTING_CHANNELS:
            msg = "Only {} channels with possibly acceptable fits. Aborting.".format(len(data))
            logger.error(msg)
            raise InversionError(msg, result=result)

    M, misfit, GFmatrix = core_inversion(data=data, **inversion_params)
    syn = GFmatrix @ full_to_deviatoric(M)

    trace_lengths = list(zip(data.index, trlen))
    result.OL2 = make_result(
        OL2Result,
        M,
        periods=(Ta, Tb),
        misfit=misfit,
        depth=event.depth,
        time_delay=t_d,
        used_traces=list(data.index),
        moment_tensor=M,
        observed_displacements=np.concatenate(data.trace),
        synthetic_displacements=syn,
        trace_lengths=OrderedDict(trace_lengths),
    )

    return data



def computeOL3(
    data: pd.DataFrame,
    event: Event,
    greens: GreensFunctions,
    result: WPhaseResult,
    processes: Optional[int] = None
):
    """Moment tensor from PDE inversion, using fixed time delay (from OL2) and
    searching for optimal centroid location + depth."""
    ol2 = result.OL2
    assert ol2 is not None
    periods = ol2.periods
    observed_displacements = ol2.observed_displacements
    t_d = ol2.time_delay
    mrf = moment_rate_function(t_d, greens.delta)
    data = data.copy()
    data["trace"] = data.apply(waveform_preprocessor(*ol2.periods), axis=1)

    logger.info("Performing grid search for best centroid location.")
    lat_grid, lon_grid = get_latlon_for_grid(
        event.latitude,
        event.longitude,
        dist_lat=3.0,
        dist_lon=3.0,
        delta=0.8
    )
    logger.debug("Grid size: %d * %d", len(lon_grid), len(lat_grid))
    latlons = [(lat, lon) for lat in lat_grid for lon in lon_grid]

    fixed_params = {
        "data": data,
        "observed_displacements": observed_displacements,
        "periods": periods,
        "t_d": t_d,
        "mrf": mrf,
        "greens": greens,
    }

    grid_search_inputs: List[dict] = [{"centroid": Point(lat, lon, event.depth)} for lat, lon in latlons]
    i_grid, _, _, grid_search_results = minimize_misfit(
        inputs=grid_search_inputs,
        parameters=fixed_params,
        processes=processes,
    )
    cenlat, cenlon = latlons[i_grid]

    logger.info("Performing grid search for best depth.")
    deps_grid = get_depths_for_grid(event.depth, greens)
    logger.debug("Depth grid size: %d", len(deps_grid))

    depth_search_inputs: List[dict] = [{"centroid": Point(cenlat, cenlon, depth)} for depth in deps_grid]
    i_dep, _, _, _ = minimize_misfit(
        inputs=depth_search_inputs,
        parameters=fixed_params,
        processes=processes
    )
    cendep = deps_grid[i_dep]

    ### Final inversion!!
    # While we already inverted at this centroid location during the depth
    # search, we threw away the corresponding synthetic waveforms; so we need
    # to re-run the inversion to get them.

    centroid = Point(latitude=cenlat, longitude=cenlon, depth=cendep)
    M, misfit, GFmatrix = core_inversion(centroid=centroid, **fixed_params)

    syn = GFmatrix @ full_to_deviatoric(M)
    trace_lengths = list(zip(data.index, data.trace.map(len)))

    result.OL3 = make_result(
        OL3Result,
        M,
        misfit=misfit,
        depth=cendep,
        time_delay=t_d,
        centroid=(cenlat, cenlon, cendep),
        periods=periods,
        used_traces=list(data.index),
        moment_tensor=M,
        observed_displacements=observed_displacements,
        synthetic_displacements=syn,
        trace_lengths=OrderedDict(trace_lengths),
        grid_search_candidates=[row["centroid"] for row in grid_search_inputs],
        grid_search_results=grid_search_results,
    )


def wpinv(
    stream: Stream,
    metadata: Mapping[str, ChannelMetadata],
    event: Event,
    gfdir: str,
    ptimes: Mapping[str, float],
    OL: int = 1,
    processes: Optional[int] = None,
    do_new_calculation: bool = False,
    ol1: Optional[OL1Result] = None,
):
    """
    This function is the main function that will compute the inversion. For
    details of the the algorithm and logic, see `here
    <https://pubs.er.usgs.gov/publication/70045140>`_.

    :param stream: The waveform data to be used as raw input
    :type stream: :py:class:`obspy.core.stream.Stream`

    :param dict metadata: Station metadata. Each key is a station ID and
        the values are instances of wphase.psi.model.ChannelMetadata.

    :param dict eqINFO: A dictionary with the preliminary event information,
        namely with keys *lat*, *lon*, *time*, *mag* (optional), and *depth*.

    :param int OL: Output level of the inversion. 1 computes just the preliminary
         magnitude, 2 perform an inversion using PDE location and 3
         uses an optimized centroid's location.

    :param OL1Result ol1:
        If provided, the preliminary magnitude calculation is skipped, with the
        provided data used instead.

    :return WPhaseResult:
    """
    # Deal with extremly shallow preliminary hypocenters by clamping to 10km
    if event.depth < 10:
        event.depth = 10

    data = build_dataframe(stream=stream, event=event, metadata=metadata, ptimes=ptimes)
    greens = GreensFunctions(gfdir)
    result = WPhaseResult(Event=event, available_traces=list(data.index))
    timer = SimpleTimer()

    if ol1 is None:
        with timer:
            computeOL1(data, result=result)
        result.timing.ol1 = timer.duration
    else:
        result.OL1 = ol1

    if OL == 1:
        return result

    with timer:
        ol2data = computeOL2(data=data, event=event, greens=greens, result=result)
    result.timing.ol2 = timer.duration

    assert result.OL2 is not None

    if OL == 2:
        return result

    if len(result.OL2.used_traces) == 0:
        logger.warning("Could not calculate OL3: no data within tolerance")
        return result

    with timer:
        computeOL3(data=ol2data, event=event, greens=greens, result=result, processes=processes)
    result.timing.ol3 = timer.duration

    if do_new_calculation:
        with timer:
            computeNew(data=data, event=event, greens=greens, result=result, t_h=result.OL2.time_delay)
        result.timing.new = timer.duration

    # for l, r, t in (
    #     ("OL1", None, result.timing.ol2),
    #     ("OL2", result.OL2.misfit, result.timing.ol2),
    #     ("OL3", result.OL3.misfit, result.timing.ol3),
    #     ("New", result.new.misfit, result.timing.new)
    # ):
    #     mf = r and f"{r:2.1f}%" or "-----"
    #     logger.info(f"{l}: {mf} in {t:3.2f}s")

    return result


def moment_rate_function(t_h, dt):
    '''
    Get the moment rate function to be convolved with the GFs.

    :param float t_h: The half duration.
    :param float dt: The time interval of the GFs.
    '''

    trianglen = 2.*int(t_h/dt)-1.
    #Five must be the minimum len for the triangle
    if trianglen <= 5:
        trianglen = 5

    mrf = triang(trianglen)
    mrf /= np.sum(mrf)

    return mrf


def get_corner_freqs_from_mag(Mw):
    '''
    Compute the corner periods required for Wphase
    inversion based on the preliminary Wphase magnitude.

    :param float Mw: Preliminary Wphase magnitude.

    :rtype: Tuple of two floats.
    '''

    if Mw >= 8.0:
        return (200., 1000.)

    # Parameters:
    Short_periods =  np.array([100, 120, 150, 200])
    Long_periods =  np.array([250, 500, 500, 1000])
    Magnitudes = np.array([6.5, 7.0, 7.5, 8.0])

    Ta_func = interp1d(Magnitudes, Short_periods, kind = 'cubic')
    Tb_func = interp1d( Magnitudes, Long_periods, kind = 'cubic')

    return (Ta_func(Mw), Tb_func(Mw))


# Attenuation relation q(Δ) for W-phase amplitudes as a function of distance.
# Lifted directly from Kanamori et al 2008.
_attenuation_relation = (
    interp1d([5.,  10., 20., 30., 40., 50.,  60.,  70.,  80.,  90.5],
             [1.5, 1.4, 1.2, 1.0, 0.7, 0.56, 0.61, 0.56, 0.50, 0.52 ],
             kind = 'cubic')
)
def correct_amplitudes_for_attenuation(amplitude: np.ndarray, distance: np.ndarray):
    return amplitude / _attenuation_relation(distance)

def preliminary_magnitude(
    tr_p2p: np.ndarray,
    dists: np.ndarray,
    azis: np.ndarray,
    trids: List[str],
    regularization: float = settings.AMPLITUDE_REGULARIZATION,
) -> PreliminaryMagnitudeFit:
    '''
    Compute the preliminary magnitude.

    :param tr_p2p: Peak to peak amplitudes for each trace.
    :param dists: station - epi distances in deg for each trace.
    :param azis: station - epi azimuths in degrees for each trace.
    :param trids: station IDs

    :return: tuple containing:

        0. Preliminary W-phase Mw magnitude.
        #. Scalar moment in dyne * cm.
        #. Prelimary estimation for the strike in degrees.
        #. Preliminary half duration for the moment rate func in seconds.
    '''

    corrected_amplitudes = correct_amplitudes_for_attenuation(tr_p2p, dists)
    azis_rad = np.array(azis) * np.pi/180. # Φ: azimuths in radians
    N = len(tr_p2p)

    # We set out to solve the system Mx = B in the least squares sense, where
    #
    # ith row of M = [1, -cos(2Φ_i), -sin(2Φ_i)]
    #            x = [2a - b, b cos(2Φ_0), b sin(2Φ_0)]
    # ith row of B = 2 max |u_w_i| / q(Δ_i)
    #
    # Applying some trig identities shows that this matrix equation is
    # equivalent to Kanamori's (12) written out for each station amplitude.
    #
    # The extra factors of 2 here are explained by our usage of p2p amplitudes
    # versus Kanamori's usage of semi-amplitudes.

    B = corrected_amplitudes
    M = np.zeros((N,3))
    M[:,0] = 1
    M[:,1] = -np.cos(2*azis_rad)
    M[:,2] = -np.sin(2*azis_rad)

    # HOWEVER, this ordinary least-squares solution is somewhat brittle: the
    # azimuthal variation of the amplitude can sometimes arrange itself such
    # that the resulting solution has b>2a and thus negative amplitudes! In
    # practice, we observed this in cases where the azimuthal coverage was very
    # poor.

    # To mitigate this, we apply regularization to penalize large values of
    # x[1] and x[2]: the cost function becomes
    #   F(x) = |Mx - B|^2 + λN(x_1^2 + x_2^2),  # N = sample count
    # which can be achieved by adding a couple rows to M and B:

    L = np.sqrt(N)*regularization
    M = np.concatenate((M, [[0, L, 0], [0, 0, L]]))
    B = np.concatenate((B, [0, 0]))

    x = lstsq(M, B, rcond=None)[0]
    amp = x[0]/2. # semi-amplitude a - b/2

    # Empirical coeff. for W-phase mag. log(a-b/2) = m*(M_w) + n
    # TODO: explain where these come from! Couldn't find them written down in
    # the paper - there was just a graph with a trendline.
    # Did someone get a ruler out!?
    m = np.log10(0.9 / 0.065)
    n = np.log10(0.9) - m * 8.5

    if amp > 0:
        # We need to add 3 because we use meters instead of mm.
        pre_wp_mag = (np.log10(amp) - n + 3)/m
    else:
        logger.warning("Preliminary magnitude fit returned a negative amplitude!")
        pre_wp_mag = 0

    unclamped = pre_wp_mag
    min_mag = settings.MINIMUM_PRELIMINARY_MAGNITUDE
    if pre_wp_mag < min_mag:
        logger.warning("Preliminary magnitude %.2f is less than %.1f. "
                       "Setting it to %.1f and continuing",
                       pre_wp_mag, min_mag, min_mag)
        pre_wp_mag = min_mag

    # Scalar moment in dyne * cm:
    M0 = 10 ** (1.5 * pre_wp_mag + 16.1)

    # Half-duration of moment rate function (just a dimensional estimate of
    # mechanism timescale):
    t_h = 1.2 * 10**-8. * M0**(1./3.)

    # Our estimate of Φ_0.
    # I need to check this. I do not know the convention
    # for the strike I should use.
    pre_strike = (0.5 * np.arctan2(x[2], x[1]) * 180./np.pi) % 360

    b = np.sqrt(x[1]*x[1] + x[2]*x[2])

    return PreliminaryMagnitudeFit(
        magnitude=pre_wp_mag,
        unclamped_magnitude=unclamped,
        scalar_moment=M0,
        timescale=t_h,
        regularization=regularization,
        strike=pre_strike,
        average_amplitude=amp,
        anisotropy=b,
        eccentricity=b/(amp+b/2),
        corrected_amplitudes=np.asarray(corrected_amplitudes),
        azimuths=np.asarray(azis),
        trids=list(trids),
    )


def fit_instrument_response(sens, freq, resp):
    """
    Fit the amplitude of the instrument response at a given set of frequencies.

    :param float sens: The instrument's sensitivity.
    :param array freq: The target frequencies.
    :param resp: The amplitudes of the response in freq.
    :type resp: array of the same length as *freq*.

    :return: Tuple containing the three time domain deconvolution coefficients.
    """

    try:
        with np.errstate(all='raise'):
            X = (2.*np.pi*freq)**2
            Y = X*X*(1./resp/resp-1.)
            fit = np.polyfit(X,Y,1)
            G = sens
            om0 = fit[1]**(.25)
            h = np.sqrt(fit[0]/4./np.sqrt(fit[1]) + 0.5)
    except FloatingPointError:
        return None
    else:
        return (om0, h, G)



def deconvolve_and_filter(
    tr: Trace,
    om0: float,
    h: float,
    G: float,
    dt: float,
    corners: int = 4,
    baselinelen: float = 60.,
    taperlen: float = 10.,
    fmin: float = 0.001,
    fmax: float = 0.005,
):
    '''
    Deconvolve and filter the trace from a single channel.

    :param tr: Data (trace) to be deconvoluted.
    :param float om0: Natural frequency of the seismometer.
    :param float h: Damping constant of the seismometer.
    :param float G: Sensitivity of the seismometer.
    :param float dt: Sampling rate in seconds.
    :param int corners: Order of the Butterworth filter to be applied.
    :param float baselinelen: Time in seconds to consider to define the
        baseline.
    :param float taperlen: Percentage of the trace to be tapered at the
        beginning.
    :param float fmin: Lower corner frequency of the Butterworth filter to be
        applied.
    :param float fmax: Upper corner frequency of the Butterworth filter to be
        applied.

    :return: Time series with the deconvoluted displacement trace
    '''

    if G <= 0.:
        # I guess this should never occur, but I was getting division by
        # zero erros in deconvolve_and_filter due either this or dt being equal to zero
        raise DeconvolutionError("Negative or zero sensitivity, skipping: {}".format(tr.id))

    if dt <= 0.:
        # I guess this should never occur, but I was getting division by
        # zero erros in deconvolve_and_filter due either this or G being equal to zero
        raise DeconvolutionError("Negative or sampling rate, skipping: {}".format(tr.id))

    if len(tr.data) <= 2:
        # Traces of length <=2 were throwing runtime exceptions; so we need to
        # explicitly skip these traces.
        raise DeconvolutionError("Trace too short to deconvolve, skipping: {}".format(tr.id))

    data = np.array(tr.data, dtype=float)

    baseline = np.mean(data[:int(baselinelen/dt)])
    nwin = int(taperlen*len(data)/100.)

    if len(data) < 2*nwin: nwin = len(data)
    data -= baseline

    taper = 0.5*(1.-np.cos(np.pi*np.linspace(0.,1.,nwin)))
    data[:nwin]*= taper

    datap1 = data[1:]
    datap2 = data[2:]

    # Acceleration + double integration
    c0 = 1./G/dt
    c1 = -2.*(1. + h*om0*dt)/G/dt
    c2 = (1. + 2.*h*om0*dt + dt*dt*om0*om0)/G/dt

    aux = c2*datap2 + c1*datap1[:-1] + c0*data[:-2]
    accel = aux[0]*np.ones(len(data))
    accel[2:] = np.cumsum(aux)

    accel = bandpassfilter(accel, dt, corners, fmin, fmax)

    vel = np.zeros(len(data))
    vel[1:] = 0.5*dt*np.cumsum(accel[:-1] + accel[1:]) # type: ignore

    dis = np.zeros(len(data))
    dis[1:] = 0.5*dt*np.cumsum(vel[:-1] + vel[1:])

    return dis


def trace_indices(traces: pd.Series, padding: int = 0) -> pd.DataFrame:
    """Given a series of traces, return a dataframe detailing the indices
    these traces will have when concatenated together."""
    trlen = traces.map(lambda tr: len(tr) + padding)
    end = np.cumsum(trlen)
    start = np.roll(end, 1)
    start[0] = 0
    return pd.DataFrame(dict(
        trlen=trlen, # this supplies the index
        start=start,
        end=end,
    ))

class Point(NamedTuple):
    latitude: float
    longitude: float
    depth: float

def get_greens_for_dataset(
    greens: GreensFunctions,
    data: pd.DataFrame,
    centroid: Point,
    periods: Tuple[float, float],
    mrf: np.ndarray,
    t_d: int,
    padding: int = 0,
):
    delta = greens.delta
    Ta, Tb = periods
    # Index of the first value that will be valid after convolution with mrf.
    t_pad = len(mrf) // 2

    # We're combining the traces of all channels into a single trace.
    # We flatten this dimension (instead of adding another numpy axis) because
    # different stations can have traces of different length.
    indices = trace_indices(data.trace, padding=padding)
    GFmatrix = np.empty((indices.end[-1], 5))

    def get_raw_synth_for(row: ChannelData):
        trlat = row.metadata.latitude
        trlon = row.metadata.longitude
        dist = locations2degrees(centroid.latitude, centroid.longitude, trlat, trlon)
        azi = gps2dist_azimuth(centroid.latitude, centroid.longitude, trlat, trlon)[1]

        # Select greens function and perform moment tensor rotation in the azimuth
        # DETAIL: The greens functions are stored for zero azimuth (which we can do
        # because of "symetry properties of a sphere". Then to recover the correct
        # moment tensor (or equivalently synthetic trace) we rotate them. For details
        # see Kanamori and Rivera 2008.
        azi_r = -azi*np.pi/180
        return greens.select_rotated(row.component, dist, centroid.depth, azi_r)

    # Stack all the raw synthetics in a single array for vectorized processing
    allsynth = np.stack(data.apply(get_raw_synth_for, axis=1)) # type: ignore

    # Use only what we think is noise to calculate and compensate the mean
    allsynth -= np.mean(allsynth[..., :60], axis=-1)[..., np.newaxis]

    # Convolve with moment rate function
    allsynth = convolve1d(allsynth, mrf, axis=-1, mode="nearest")[..., t_pad-2:-t_pad]

    # Pad data back to original length by extending values at edges
    allsynth = pad1d(allsynth, (t_pad, t_pad), axis=-1, mode="edge")
    allsynth -= np.mean(allsynth[..., :60], axis=-1)[..., np.newaxis]

    # Bandpass filter to extract W phase
    allsynth = bandpassfilter(allsynth, delta, 4, 1/Tb, 1/Ta, axis=-1)

    for synth, (_, row) in zip(allsynth, (data.join(indices)).iterrows()):
        # The time interval we copy here depends on each trace's p-time, thus the loop.
        synth = ltrim(synth, row.ptime - t_d, delta)
        synth = synth[:, :row.trlen]

        output = GFmatrix[row.start:row.end]

        # the first of the following lines are due to constraint that volume does not change
        output[:, 0] = synth[0][:] - synth[1][:]
        output[:, 1] = synth[2][:] - synth[1][:]
        output[:, 2] = synth[4][:]
        output[:, 3] = synth[5][:]
        output[:, 4] = synth[3][:]

    return GFmatrix

class InversionResult(NamedTuple):
    tensor: np.ndarray
    misfit: float
    greens_matrix: np.ndarray

def _inversion_for_scipy(
    x: Sequence,
    data: pd.DataFrame,
    periods: Tuple[float, float],
    greens: GreensFunctions,
    t_d: float,
    t_h: float,
    observed_displacements: Optional[np.ndarray] = None,
) -> InversionResult:
    """core_inversion wrapper for use with scipy.optimize."""
    lat, lon, depth = x
    return core_inversion(
        t_d=round(t_d),
        centroid=Point(lat, lon, depth),
        data=data,
        periods=periods,
        mrf=moment_rate_function(t_h, greens.delta),
        greens=greens,
        observed_displacements=observed_displacements,
    )

def _misfit_for_scipy(
    x: Sequence,
    data: pd.DataFrame,
    periods: Tuple[float, float],
    greens: GreensFunctions,
    t_d: float,
    t_h: float,
    observed_displacements: Optional[np.ndarray] = None,
) -> float:
    return _inversion_for_scipy(x, data, periods, greens, t_d, t_h, observed_displacements).misfit

class NonlinearParams(NamedTuple):
    t_d: float
    t_h: float
    latitude: float
    longitude: float
    depth: float

def omni_minimize(
    data: pd.DataFrame,
    greens: GreensFunctions,
    periods: Tuple[float, float],
    t_d: float,
    t_h: float,
    x0: Point,
    tol: Optional[float] = 1e-2,
) -> Tuple[Point, InversionResult]:
    if "trace" not in data:
        data = data.copy()
        data["trace"] = data.apply(waveform_preprocessor(*periods), axis=1)
    misfit = _misfit_for_scipy(x0, data, periods, greens, t_d, t_h)
    logger.info("Initial guess: (%.2f %.2f %.2f) => %.1f%%", *x0, misfit)
    obs = np.concatenate(data.trace)
    optimized = minimize(
        _misfit_for_scipy,
        x0,
        args=(
            data,
            periods,
            greens,
            t_d,
            t_h,
            obs,
        ),
        method="Nelder-Mead",
        bounds=(
            (-90, 90),
            (None, None),
            (greens.depths.min(), greens.depths.max()),
        ),
        options={
            "xatol": tol,
            "fatol": tol,
            "initial_simplex": [
                [x0[0], x0[1], 4*x0[2]],
                [x0[0]-4, x0[1]-2, 1],
                [x0[0]+4, x0[1]-2, 1],
                [x0[0], x0[1]+4, 1],
            ],
        },
    )
    inversion = _inversion_for_scipy(optimized.x, data, periods, greens, t_h, t_h, obs)
    logger.info("Minimizer: (%.2f %.2f %.2f) => %.1f%%", *optimized.x, optimized.fun)
    return Point(*optimized.x), inversion

def core_inversion(
    data: pd.DataFrame,
    periods: Tuple[float, float],
    t_d: int,
    centroid: Point,
    mrf: np.ndarray,
    greens: GreensFunctions,
    observed_displacements: Optional[np.ndarray] = None,
) -> InversionResult:
    '''
    Perform the actual W-phase inversion.

    :param array observed_displacements: Array containing concatenated traces of observed displacements. If omitted, we rebuild it from data.
    :param float t_d: Time delay
    :param tuple cmtloc: (lat, lon, depth) of the centroid's location.
    :param tuple periods: (Ta,Tb), passband periods.
    :param array mrf: Moment rate function.
    :param array trlen: Array containing the length of each trace.
    :param dict metadata: Dictionary with the metadata of each station.
    :param list trlist: List with the station id which will contribute to the inv.
    :param greens: Path to the greens functions or a GreensFunctions instances.
    '''

    GFmatrix = get_greens_for_dataset(
        greens=greens,
        data=data,
        centroid=centroid,
        periods=periods,
        mrf=mrf,
        t_d=t_d,
    )

    # perform the inversion
    obs: np.ndarray = (
        np.concatenate(data.trace)
        if observed_displacements is None
        else observed_displacements
    )
    M, residual, *_  = lstsq(GFmatrix, obs, rcond=None)

    # def _log(label, x):
    #     logger.info(f"{label}: {x.dtype}{x.shape}\n{x.flags}")
    #_log("GF", GFmatrix)
    #_log("observed", observed_displacements)
    #_log("M", M)

    # construct the synthetics
    syn = GFmatrix @ M

    # from numpy.testing import assert_almost_equal
    # assert_almost_equal(np.sum((syn-observed_displacements)**2), residual)
    # assert_almost_equal(np.sum(syn**2), syn2)
    #syn = GFmatrix @ M
    misfit = 100 * np.sqrt(residual / (syn @ syn))
    # misfit = compute_misfit(residual, GFmatrix, M)

    logger.debug(
        "core_inversion(t_d=%d, c=(%.2f, %.2f, %.2f), #tr=%d) => %.1f%%",
        t_d, centroid.latitude, centroid.longitude, centroid.depth, len(data), misfit
    )

    return InversionResult(tensor=deviatoric_to_full(M), misfit=misfit, greens_matrix=GFmatrix)

def deviatoric_to_full(mt: np.ndarray) -> np.ndarray:
    """Convert the 5-element representation of a deviatoric moment tensor to
    its 6-element representation."""
    return np.insert(mt, 2, - mt[0] - mt[1])

def full_to_deviatoric(mt: np.ndarray) -> np.ndarray:
    """Convert the 6-element representation of a deviatoric moment tensor to
    its 5-element representation."""
    return np.delete(mt, 2)

def compute_misfit(residual, GFmatrix, M):
    syn2 = np.einsum("ij,j,ik,k->", GFmatrix, M, GFmatrix, M,
                     optimize=['einsum_path', (0, 1), (0, 2), (0, 1)])
    return 100*np.sqrt(residual/syn2)


def get_depths_for_grid(hypdep, greens, length=100):
    '''
    This function will return a list with the adequate values
    of depth to look at in the grid search
    '''

    ####Parameters:
    min_dep = 11. #No shallower depths allowed.
     # Total length of the section in which we'll search
    #####
    depths = np.asarray(greens.depths, dtype=float)
    depths_grid = depths[np.where(np.abs(depths-hypdep)<length*0.5)]
    depths_grid = depths_grid[np.where(depths_grid > min_dep)]

    return depths_grid



def get_latlon_for_grid(hyplat,hyplon, dist_lat = 1.2,
                            dist_lon = 1.2,delta = 0.4 ):
    '''
    This function will return a list with the adequate values
    of lan and lon to look at in the grid search

    :param float dist_lat: Half distance in degrees to search.
    :param float dist_lon: Half distance in degrees to search.
    :param float delta: Spacing of the grid in degrees.

    :return: Tuple containing the latitudes and longitudes to
        include in the search as numpy arrays.
    '''

    lat_grid = np.arange(hyplat-dist_lat,hyplat+dist_lat,delta)
    lon_grid = np.arange(hyplon-dist_lon,hyplon+dist_lon,delta)

    return lat_grid, lon_grid



def core_inversion_wrapper(kwargs: Mapping[str, Any], fixed_kwargs: Mapping[str, Any]):
    '''
    Wrapper for :py:func:`core_inversion` for use with
    :py:class:`multiprocessing.Pool`.

    Doesn't return greens functions to avoid pickling overhead.
    '''
    result = core_inversion(**{**fixed_kwargs, **kwargs})
    return result.tensor, result.misfit


def minimize_misfit(inputs: Sequence[Mapping], parameters: Mapping[str, Any], processes=None):
    """Given a list of input tuples for core_inversion, find the one with the
    lowest misfit.

    :param inputs: an iterable of mappings to be provided to core_inversion as kwargs.
    :param parameters: fixed kwargs for core_inversion
    :param processes: number of parallel processes to use.
    :return: A tuple (index_of_minimizer, moment_tensor, misfit, results list)."""
    worker = partial(core_inversion_wrapper, fixed_kwargs=parameters)
    # with ProcessPoolExecutor(max_workers=processes) as pool:
    #     results = list(pool.map(worker, inputs))
    results = list(map(worker, inputs))
    misfits = np.array([x[1] for x in results])
    i_min = misfits.argmin()
    return i_min, results[i_min][0], results[i_min][1], results


def get_timedelay_misfit(t_d, GFmatrix, trlen, observed_displacements, max_t_d):
    max_t_d = int(max_t_d)
    t_d = int(t_d)
    t_d2 = max_t_d - t_d
    GFmatrix_sm = np.zeros((np.array(trlen, dtype=int).sum(),5))
    cumtrlens = np.concatenate(([0], np.array(trlen, dtype=int).cumsum()))
    tb_fm = 0
    for i_ntr, ta in enumerate(cumtrlens[:-1]):
        l = trlen[i_ntr]
        tb = ta + l
        ta_fm = tb_fm
        tb_fm = ta_fm +  max_t_d + tb-ta
        GFmatrix_sm[ta:tb] = GFmatrix[ta_fm + t_d2 :ta_fm + t_d2 + l]

    inversion = lstsq(GFmatrix_sm, observed_displacements, rcond=None)
    return inversion[1][0]


def make_result(cls, M: ArrayLike, misfit: float, **extra):
    """Given a moment tensor, print some pretty log messages describing it and
    convert it to a standard result object.

    Parameters
    ----------
    cls :
        Result class to construct: OL2Result or OL3Result
    M :
        Moment tensor as a 6-element list
    misfit :
        Inversion misfit as a percentage
    extra :
        Extra args to pass to the result constructor

    Return
    ------
    OL2Result
    """
    M = np.asarray(M)
    logger.info("%s Result:" % cls.__name__)
    logger.info("%5s % 10s % 10s % 10s", "", "r", "t", "p")
    logger.info("%5s %+.3e %+.3e %+.3e", "r", M[0], M[3], M[4])
    logger.info("%5s % 10s %+.3e %+.3e", "t", "",   M[1], M[5])
    logger.info("%5s % 10s % 10s %+.3e", "p", "",   "",   M[2])
    logger.info("misfit: %.0f%%", misfit)
    M2 = M**2
    m0 = np.sqrt(0.5 * (M2[0] + M2[1] + M2[2]) + M2[3] + M2[4] + M2[5])
    mag = 2./3.*(np.log10(m0)-9.10)
    logger.info("m0: % e", m0)
    logger.info("magnitude: %.5f", mag)
    return cls(
        Mrr=M[0],
        Mtt=M[1],
        Mpp=M[2],
        Mrt=M[3],
        Mrp=M[4],
        Mtp=M[5],
        misfit=misfit,
        m0=m0,
        magnitude=round(mag, 1),
        **extra,
    )
