# -*- coding: utf-8 -*-
"""Core implementation of the W-Phase inversion."""
from __future__ import absolute_import, print_function

from builtins import str
from builtins import range
from collections import OrderedDict
import sys, os, glob, traceback, logging
from concurrent.futures import ProcessPoolExecutor
from timeit import default_timer as timer
from typing import Sequence
import numpy as np
from numpy.linalg import lstsq
from scipy.interpolate import interp1d
from scipy.signal import triang
from scipy import ndimage
from scipy.signal import detrend
import h5py
from obspy.core import read, Stream, UTCDateTime

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
from .greens import GreensFunctions
from .bandpass import bandpassfilter
from .model import (
    WPhaseResult,
    OL1Result,
    OL2Result,
    OL3Result,
    Event,
    TimeDelayMisfits,
)
from .exceptions import InversionError, RTdeconvError
from .seismoutils import azimuthal_gap, ltrim, rot_12_NE

logger = logging.getLogger(__name__)

def wpinv(
    st: Stream,
    metadata: dict,
    event: Event,
    gfdir: str,
    OL: int = 1,
    processes: int = None,
):
    """
    This function is the main function that will compute the inversion. For
    details of the the algorithm and logic, see `here
    <https://pubs.er.usgs.gov/publication/70045140>`_.

    :param st: The waveform data to be used as raw input
    :type st: :py:class:`obspy.core.stream.Stream`

    :param dict metadata: Station metadata. Each key is a station ID and
        the values are a dictionaries with the metadata. Each dictionary
        should look like:

        .. code-block:: python

            {'azimuth': 0.0,
                'dip': -90.0,
                'elevation': 19.0,
                'gain': 59680600.0,
                'latitude': 2.0448,
                'longitude': -157.4457,
                'poles': [
                    (-0.035647-0.036879j),
                    (-0.035647+0.036879j),
                    (-251.33+0j),
                    (-131.04-467.29j),
                    (-131.04+467.29j)],
                'sensitivity': 33554000000.0,
                'transfer_function': 'A',
                'zeros': [0j, 0j]}

    :param dict eqINFO: A dictionary with the preliminary event information,
        namely with keys *lat*, *lon*, *time*, *mag* (optional), and *depth*.

    :param int OL: Output level of the inversion. 1 computes just the preliminary
         magnitude, 2 perform an inversion using PDE location and 3
         uses an optimized centroid's location.

    :return: What is returned depends on the processing level (*OL*), as described below.

        - OL = 1, a tuple with elements:
            0. Preliminary Mw magnitude of the event.
            #. List with the stations contributing to the preliminary
                magnitude. (epicentral order)
            #. List with the peak-to-peak amplitude (meters) of each station
                contributing to the preliminary magnitude (epicentral order).

        - OL = 2, a tuple with elements:
            0. Six component of the moment tensor in Nm,
                ['RR', 'PP', 'TT', 'TP' , 'RT', 'RP'] (ready to be plotted with
                Beachball from ObsPy).
            #. Array with the concatenated traces of observed displacements
                sorted according to epicentral distance.
            #. Array with the concatenated traces of synthetic seismograms
                sorted according to epicentral distance.
            #. List with the stations contributing to the final solution,
                sorted according to epicentral distance.
            #. List with the lengths of each trace in trlist. Note that
                sum(trlen) = len(observed_displacements) = len(syn).

        - OL = 3, Same as OL2 plus:
            5. Optimal centroid location (lat, lon, dep).
            #. Time delay/Half duration of the source (secs).
            #. latitude-longitude grid used to find the optimal centroid
                location
            #. Inversion result for each grid point in latlons
            #. Dictionary with relavant information about the data processing
                so it can redone easily outside this function.
    """
    #############Output Level 1#############################
    #############Preliminary magnitude:#####################

    # Ta and Tb give the corner periords in seconds of the band pass filter
    # used at this stage of processing.
    Ta = 200.
    Tb = 1000.

    # vector of periods to consider in the instrument response fitting.
    T = np.linspace(Ta,Tb,500)
    freq = 1./T

    # convert to angular frequency
    omega = freq*2.*np.pi

    # create a function to compute the amplitude of the response from poles and
    # zeros.
    Vpaz2freq = np.vectorize(paz_2_amplitude_value_of_freq_resp)

    hyplat = event.latitude
    hyplon = event.longitude
    hypdep = event.depth
    orig = UTCDateTime(event.time)

    # Deal with extremly shallow preliminary hypocenter
    if hypdep < 10.:
        hypdep = 10.

    #In case of multiple locations we favor
    # '00'. This may be improved if we select the one with the longer period
    # sensor.
    st_sel = st.select(location = '00')
    st_sel+= st.select(location = '--')
    st_sel+= st.select(location = '') #Check this.

    # We don't really care exactly what ordering obspy chooses
    # here; we just want a deterministic ordering so we can compare runs.
    st_sel.sort()

    #st_sel = st_sel.select(component = 'Z')
    #We also want to select stations with one location code which is not the
    #default (see above).

    logger.info('Initial number of traces: using %d with default locations (of %d total)', len(st_sel), len(st))

    # rotate the horizontal components to geographical north or east
    st_sel = rot_12_NE(st_sel, metadata)

    # traces to use the preliminary mag calculation
    st_sel_prem_mag = st_sel.select(component = 'Z').copy()

    # will contain distances epicenter - station
    DIST = np.array([])

    # trlist  contains the station IDs. If a station is then rejected
    # it must be removed from this list.
    trlist = [tr.id for tr in st_sel]
    trlist_pre = []

    # Peak to peak amplitude of each trace.
    tr_p2p = []

    # List with the station azimuths
    AZI = []

    # Instrument deconvolution and bandpass filtering:
    i = 0
    for tr in st_sel_prem_mag:
        trmeta =  metadata[tr.id]
        trlat = trmeta['latitude']
        trlon = trmeta['longitude']
        tr.data = np.array(tr.data, dtype=float)

        # compute distance in degrees between 2 locations
        dist = locations2degrees(hyplat, hyplon, trlat, trlon)

        # compute azimuth from north
        azi = gps2dist_azimuth(hyplat, hyplon, trlat, trlon)[1]

        # if p-arrival time is not calculated, then calculate it.
        # WARNING: this can be very slow in new versions of obspy.
        t_p = trmeta.get('ptime')
        if not t_p:
            from obspy.taup.taup import getTravelTimes
            t_p =  getTravelTimes(dist,hypdep)[0]['time']

        # sample period in seconds
        dt = tr.stats.delta

        # Wphase (UTC) time window
        t1 = orig + t_p
        t2 = t1 + dist*settings.WPHASE_CUTOFF

        # accounting for the units of the transfer function in the instrument response... see README.
        if trmeta['transfer_function'] == "B":
            AmpfromPAZ  = Vpaz2freq(trmeta,freq/2./np.pi)  # hz to rad*hz
        elif trmeta['transfer_function'] == "A":
            AmpfromPAZ  = Vpaz2freq(trmeta,freq)
        else:
            logger.warning("Unknown transfer function. Skipping %s", tr.id)
            trlist.remove(tr.id)
            continue

        # Fitting the instrument response and getting coefficients
        response_coefficients = fit_instrument_response(
            trmeta['sensitivity'], freq, AmpfromPAZ
        )
        if response_coefficients is None:
            logger.warning("Impossible to get Coeff. Skipping %s", tr.id)
            trlist.remove(tr.id)
            continue
        else:
            om0, h, G = response_coefficients

        AmpfromCOEFF= np.abs(omega*omega / \
                (omega*omega + 2j*h*om0*omega - om0*om0))

        # L2 norm:
        misfit = 100*np.linalg.norm(AmpfromPAZ-AmpfromCOEFF) \
                / np.linalg.norm(AmpfromPAZ)

        if misfit > settings.RESPONSE_MISFIT_TOL:
            logger.warning('Bad fitting for response (misfit=%e). Skipping %s', misfit, tr.id)
            continue

        # tr.data will cointain the deconvolved and filtered displacements.
        try:
            tr.data, coeff = RTdeconv(
                tr,
                om0,
                h,
                G,
                dt,
                corners=4,
                baselinelen=60./dt,
                taperlen=10.,
                fmin=1./Tb,
                fmax=1./Ta,
                get_coef=True)

        except RTdeconvError as e:
            logger.warning("Error deconvolving trace %s: %s", tr.id, str(e))
            trlist.remove(tr.id)
            continue

        # trim to the Wphase time window
        tr.trim(t1,t2)
        if len(tr)== 0:
            logger.warning("Empty trace. Skipping %s", tr.id)
            trlist.remove(tr.id)
            continue

        trlist_pre.append(tr.id)
        tr_p2p.append(tr[:].max()-tr[:].min())
        AZI.append(azi)
        DIST = np.append(DIST, dist)
        i += 1

    # Sorting the IDs according to their distance to the source:
    sorted_indices = np.argsort(DIST)
    trlist_pre = [trlist_pre[i] for i in sorted_indices]
    tr_p2p = [tr_p2p[i] for i in sorted_indices]
    AZI = [AZI[i] for i in sorted_indices]
    DIST = np.sort(DIST)

    # Median rejection
    median_p2p = np.nanmedian(tr_p2p)
    mrcoeff = settings.MEDIAN_REJECTION_COEFF
    p2pmax, p2pmin = mrcoeff[1]*median_p2p, mrcoeff[0]*median_p2p
    accepted_traces = []
    for i, (trid, p2p) in enumerate(zip(trlist_pre, tr_p2p)):
        if np.isnan(p2p):
            logger.warning("P2P Amp for %s is NaN! Excluding.", trid)
        elif p2p > p2pmax:
            logger.info("P2P Amp for %s is too big (%.2e > %.2e). Excluding.", trid, p2p, p2pmax)
        elif p2p < p2pmin:
            logger.info("P2P Amp for %s is too small (%.2e < %.2e). Excluding.", trid, p2p, p2pmin)
        else:
            accepted_traces.append(i)

    gap = azimuthal_gap(AZI)
    if gap > settings.MAXIMUM_AZIMUTHAL_GAP:
        raise InversionError("Lack of azimuthal coverage (%.0f° > %.0f°). Aborting."
                             % (gap, settings.settings.MAXIMUM_AZIMUTHAL_GAP))
    if len(accepted_traces) < settings.MINIMUM_STATIONS:
        raise InversionError("Lack of stations (%d < %d). Aborting."
                             % (len(accepted_traces), settings.MINIMUM_STATIONS))

    logger.info("Traces accepted for preliminary magnitude calculation: {}"
                .format(len(accepted_traces)))

    # Get the preliminary mag:
    tr_p2p_con = [tr_p2p[i] for i in accepted_traces]
    DIST_con   = [DIST[i] for i in accepted_traces]
    AZI_con    = [AZI[i] for i in accepted_traces]
    trlist_pre_con = [trlist_pre[i] for i in accepted_traces]
    pre_results = preliminary_magnitude(tr_p2p_con, DIST_con, AZI_con, trlist_pre_con)

    pre_wp_mag = pre_results['magnitude']
    t_h = pre_results['t_h']
    logger.info("Preliminary t_h = %.2f" % t_h)

    logger.info("OL1:")
    logger.info("Average amplitude:  %.1em", pre_results['average_amplitude'])
    logger.info("Magnitude:          %.2f", pre_wp_mag)
    logger.info("Strike:             %.1f°", pre_results['strike'])
    logger.info("Eccentricity (b/a): %.2f", pre_results['eccentricity'])

    ol1 = OL1Result(
        preliminary_calc_details=pre_results,
        used_traces=trlist,
        nstations=len(accepted_traces),
        magnitude=round(pre_wp_mag, 1),
    )
    result = WPhaseResult(OL1=ol1, Event=event)
    if OL == 1:
        return result
    #############  Output Level 2    #######################################
    #############  Moment Tensor based on preliminary hypercenter (PDE) ####
    # Much of what follows is the same as what was done above, but we will be
    # using a different set of stations and can handle multiple components per
    # station (i.e. horizontals also)
    ########################################################################

    # Redefine and define some values according to the pre_wp_mag
    Ta, Tb = get_corner_freqs_from_mag(pre_wp_mag)
    logger.info("Filter corner periods: %.1f, %1f" % (Ta, Tb))
    T = np.linspace(Ta,Tb,500)
    freq = 1./T
    omega = freq*2.*np.pi

    # Build the Moment Rate Function (MRF):
    greens = GreensFunctions(gfdir)
    dt = greens.delta
    MRF = MomentRateFunction(t_h, dt)

    # Building a stream with the synthetics and the observed displacements vector
    tr_p2p = [] #Peak to peak amplitude of each trace.

    # Preparing the data:
    DIST = []
    DATA_INFO = {} #Minimum info to be able to filter the displacements afterwards
    for itr, trid in enumerate(trlist[:]):
        trmeta =  metadata[trid]
        trlat = trmeta['latitude']
        trlon = trmeta['longitude']
        dist = locations2degrees(hyplat, hyplon, trlat, trlon)
        t_p = trmeta.get('ptime')
        if not t_p:
            from obspy.taup.taup import getTravelTimes
            t_p =  getTravelTimes(dist,hypdep)[0]['time']
        dt = tr.stats.delta

        # W-phase time window
        t1 = orig + t_p
        t2 = t1 + dist*settings.WPHASE_CUTOFF
        tr = st_sel.select(id = trid)[0]

        tr.data = np.array(tr.data, dtype=float)
        tf = trmeta['transfer_function']
        if tf == "B":
            AmpfromPAZ  = Vpaz2freq(trmeta,freq/2./np.pi)  # hz to rad*hz
        elif tf == "A":
            AmpfromPAZ  = Vpaz2freq(trmeta,freq)
        else:
            logger.warning("%s: Unknown transfer function %s. Skipping this trace", tr.id, tf)
            trlist.remove(tr.id)
            continue

        response_coefficients = fit_instrument_response(
            trmeta['sensitivity'], freq, AmpfromPAZ
        )
        if response_coefficients is None:
            logger.warning("%s: Could not fit instrument response. Skipping this trace", tr.id)
            trlist.remove(tr.id)
            continue
        else:
            om0, h, G = response_coefficients

        DATA_INFO[tr.id] = [om0, h, G, dt, t1, t2]

        try:
            tr.data, coeff = RTdeconv(
                tr, om0, h, G, dt,
                corners=4,
                baselinelen=60./dt,
                taperlen= 10.,
                fmin = 1./Tb,
                fmax = 1./Ta,
                get_coef = True)

        except RTdeconvError as e:
            logger.warning("%s: Skipping due to error in deconvolution: %s", tr.id, str(e))
            trlist.remove(tr.id)
            continue

        #Check length of the trace:
        tr.trim(t1, t2)
        if len(tr) == 0:
            logger.warning("%s: Trace is empty. Skipping this trace", tr.id)
            trlist.remove(tr.id)
            continue
        tr_p2p.append(tr[:].max() - tr[:].min())

        DIST.append(dist)

    DIST = np.array(DIST)
    trlist = [trlist[i] for i in np.argsort(DIST)]
    tr_p2p = [tr_p2p[i] for i in np.argsort(DIST)]

    # Rejecting outliers:
    observed_displacements = np.array([]) # observed displacements vector.
    trlen = [] # A list with the length of each station data.
    trlist2 = trlist[:]

    mrcoeff = settings.MEDIAN_REJECTION_COEFF
    median_AMP = np.median(tr_p2p)
    for i, amp in enumerate(tr_p2p):
        if (amp > median_AMP*mrcoeff[1]
        or amp < median_AMP*mrcoeff[0]):
            trlist2.remove(trlist[i])
            logger.warning("%s: Amplitude is %.2fx the median, which is "
                           "outside our acceptable range of [%.2f, %.2f]. "
                           "Rejecting this trace.",
                           trlist[i],
                           amp/median_AMP,
                           mrcoeff[0],
                           mrcoeff[1])

        else:
            tr = st_sel.select(id = trlist[i])[0]
            observed_displacements =  np.append(observed_displacements, tr.data[:],0)
            trlen.append(len(tr))
    trlist = trlist2[:]

    logger.info('number of traces for OL2: %d', len(trlist))

    #### Inversion:
    # Search for the optimal t_d
    ##I need to test the optimal way to do this. map is taking too much time

    # time delays for which we will run inversions (finally choosing the one with
    # lowest misfit)
    time_delays = np.arange(1., settings.MAXIMUM_TIME_DELAY)

    # extract the greens function matrix for all time delays. Note that this does not
    # perform an inversion, because OnlyGetFullGF=True.
    GFmatrix = core_inversion(
        0,
        (hyplat, hyplon, hypdep),
        (Ta, Tb),
        MRF,
        observed_displacements,
        trlen,
        metadata,
        trlist,
        gfdir=gfdir,
        OnlyGetFullGF=True,
        max_t_d=settings.MAXIMUM_TIME_DELAY,
    )

    # inputs for the TIME DELAY search
    inputs = [(t_d_test, GFmatrix, trlen, observed_displacements, settings.MAXIMUM_TIME_DELAY) for t_d_test in time_delays]
    with ProcessPoolExecutor(max_workers=processes) as pool:
        misfits = list(pool.map(get_timedelay_misfit_wrapper,inputs))

    # Set t_d (time delay) and and t_h (half duration) to optimal values:
    mis_min = int(np.array(misfits).argmin())
    result.misfits = TimeDelayMisfits(array=misfits, min=mis_min)
    t_d = t_h = time_delays[mis_min]
    MRF = MomentRateFunction(t_h, dt)
    logger.info("Source time function, time delay: %d, %f", len(MRF), t_d)
    logger.info("revised t_h = %.2f" % t_h)

    #### Removing individual bad fitting. this recursively removes stations with misfits
    # outside of the acceptable range defined by the setting MISFIT_TOL_SEQUENCE.
    M, misfit, GFmatrix = core_inversion(
        t_d,
        (hyplat, hyplon, hypdep),
        (Ta, Tb),
        MRF,
        observed_displacements,
        trlen,
        metadata,
        trlist,
        gfdir=gfdir,
        return_gfs=True)

    # Remove bad traces
    for tol in settings.MISFIT_TOL_SEQUENCE:
        # Make sure there are enough channels
        GFmatrix, observed_displacements, trlist, trlen = remove_individual_traces(
            tol,
            M,
            GFmatrix,
            observed_displacements,
            trlist,
            trlen)

        if len(trlist) < settings.MINIMUM_FITTING_CHANNELS:
            msg = "Only {} channels with possibly acceptable fits. Aborting.".format(len(trlist))
            logger.error(msg)
            raise InversionError(msg)

        M, misfit, GFmatrix = core_inversion(
            t_d,
            (hyplat, hyplon, hypdep),
            (Ta,Tb),
            MRF,
            observed_displacements,
            trlen,
            metadata,
            trlist,
            gfdir=gfdir,
            return_gfs=True
        )

    syn = (M[0]*GFmatrix[:,0] + M[1]*GFmatrix[:,1] +
           M[3]*GFmatrix[:,2] + M[4]*GFmatrix[:,3] +
           M[5]*GFmatrix[:,4])


    trace_lengths = list(zip(trlist, trlen))
    result.OL2 = make_result(
        OL2Result,
        M,
        misfit=misfit,
        depth=hypdep,
        time_delay=t_d,
        used_traces=trlist,
        moment_tensor=M,
        observed_displacements=observed_displacements,
        synthetic_displacements=syn,
        trace_lengths=OrderedDict(trace_lengths),
    )

    if len(trlen) == 0:
        logger.warning("Could not calculate OL3: no data within tolerance")
        return result
    elif OL == 2:
        return result

    #############  Output Level 3   #############################
    ###Moment Tensor based on grid search  ######################
    ###for the optimal centroid's location #####################

    logger.info("Performing grid search for best centroid location.")
    lat_grid, lon_grid = get_latlon_for_grid(hyplat, hyplon, dist_lat=3.0,
                                             dist_lon=3.0, delta=0.8)
    logger.debug("Grid size: %d * %d", len(lon_grid), len(lat_grid))
    latlons = [(lat, lon) for lat in lat_grid for lon in lon_grid]

    grid_search_inputs = [
        (t_d, (lat, lon, hypdep), (Ta, Tb), MRF,
         observed_displacements, trlen, metadata, trlist, greens)
        for lat, lon in latlons
    ]

    i_grid, _, _, grid_search_results = minimize_misfit(
        grid_search_inputs,
        processes=processes
    )
    cenlat, cenlon = latlons[i_grid]

    logger.info("Performing grid search for best depth.")
    deps_grid = get_depths_for_grid(hypdep, greens)
    logger.debug("Depth grid size: %d", len(deps_grid))

    depth_search_inputs = [
        (t_d, (cenlat, cenlon, depth), (Ta,Tb),
         MRF, observed_displacements, trlen, metadata, trlist, greens)
        for depth in deps_grid
    ]

    i_dep, _, _, _ = minimize_misfit(
        depth_search_inputs,
        processes=processes
    )
    cendep = deps_grid[i_dep]

    ### Final inversion!!
    # While we already inverted at this centroid location during the depth
    # search, we threw away the corresponding synthetic waveforms; so we need
    # to re-run the inversion to get them.

    M, misfit, GFmatrix = core_inversion(t_d, (cenlat, cenlon, cendep),
                                         (Ta,Tb), MRF, observed_displacements,
                                         trlen, metadata, trlist, gfdir=gfdir,
                                         return_gfs=True)

    syn = (M[0]*GFmatrix[:,0] + M[1]*GFmatrix[:,1]
          + M[3]*GFmatrix[:,2] + M[4]*GFmatrix[:,3]
          + M[5]*GFmatrix[:,4])


    trace_lengths = list(zip(trlist, trlen))

    result.OL3 = make_result(
        OL3Result,
        M,
        misfit=misfit,
        depth=cendep,
        time_delay=t_d,
        centroid=(cenlat, cenlon, cendep),
        used_traces=trlist,
        moment_tensor=M,
        observed_displacements=observed_displacements,
        synthetic_displacements=syn,
        trace_lengths=OrderedDict(trace_lengths),
        grid_search_candidates=[row[1] for row in grid_search_inputs],
        grid_search_results=grid_search_results,
    )

    result.CreationTime = UTCDateTime.now()
    return result


def MomentRateFunction(t_h, dt):
    '''
    Get the moment rate function to be convolved with the GFs.

    :param float t_h: The half duration.
    :param float dt: The time interval of the GFs.
    '''

    trianglen = 2.*int(t_h/dt)-1.
    #Five must be the minimum len for the triangle
    if trianglen <= 5:
        trianglen = 5

    MRF = triang(trianglen)
    MRF /= np.sum(MRF)

    return MRF


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



def preliminary_magnitude(tr_p2p, dists, azis, trids,
                          regularization=settings.AMPLITUDE_REGULARIZATION):
    '''
    Compute the preliminary magnitude.

    :param list tr_p2p: Peak to peak amplitudes for each trace.
    :param list dists: station - epi distances in deg for each trace.
    :param list azis: station - epi azimuths in degrees for each trace.
    :param list trids: station IDs

    :return: tuple containing:

        0. Preliminary W-phase Mw magnitude.
        #. Scalar moment in dyne * cm.
        #. Prelimary estimation for the strike in degrees.
        #. Preliminary half duration for the moment rate func in seconds.
    '''

    # Attenuation relation q(Δ) for W-phase amplitudes as a function of distance.
    # Lifted directly from Kanamori et al 2008.
    q = interp1d([5.,  10., 20., 30., 40., 50.,  60.,  70.,  80.,  90.5],
                 [1.5, 1.4, 1.2, 1.0, 0.7, 0.56, 0.61, 0.56, 0.50, 0.52 ],
                 kind = 'cubic')
    corrected_amplitudes = tr_p2p / q(dists)

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
        pre_wp_mag = 6.5

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

    return dict(
        magnitude=pre_wp_mag,
        unclamped_magnitude=unclamped,
        M0=M0,
        t_h=t_h,
        regularization=regularization,
        strike=pre_strike,
        average_amplitude=amp,
        anisotropy=b,
        eccentricity=b/(amp+b/2),
        corrected_amplitudes=corrected_amplitudes,
        azimuths=azis,
        trids=trids,
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



def RTdeconv(
        tr,
        om0,
        h,
        G,
        dt,
        corners = 4,
        baselinelen = 60.,
        taperlen = 10.,
        fmin = 0.001,
        fmax = 0.005,
        get_coef = False,
        data_type = 'raw'):
    '''
    Using the coef. of one station computes the displacement trace.

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
    :param str data_type: Type of data. 'raw' for raw data to be deconvoluted,
        'accel' for deconvoluted acceleration.

    :return: Time series with the deconvoluted displacement trace
    '''

    if G <= 0.:
        # I guess this should never occur, but I was getting division by
        # zero erros in RTdeconv due either this or dt being equal to zero
        raise RTdeconvError("Negative or zero sensitivity, skipping: {}".format(tr.id))

    if dt <= 0.:
        # I guess this should never occur, but I was getting division by
        # zero erros in RTdeconv due either this or G being equal to zero
        raise RTdeconvError("Negative or sampling rate, skipping: {}".format(tr.id))

    if len(tr.data) <= 2:
        # Traces of length <=2 were throwing runtime exceptions; so we need to
        # explicitly skip these traces.
        raise RTdeconvError("Trace too short to deconvolve, skipping: {}".format(tr.id))

    data = np.array(tr.data,dtype='float')
    sta  = tr.stats.station

    if data_type == 'raw':
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
    elif data_type == 'accel':
        accel = data

    accel = bandpassfilter(accel, dt, corners, fmin, fmax)

    vel = np.zeros(len(data))
    vel[1:] = 0.5*dt*np.cumsum(accel[:-1]+accel[1:])

    dis = np.zeros(len(data))
    dis[1:] = 0.5*dt*np.cumsum(vel[:-1]+vel[1:])

    if get_coef:
        return dis, (c0,c1,c2)

    return dis


def core_inversion(
    t_d,
    cmtloc,
    periods,
    MRF,
    observed_displacements,
    trlen,
    metadata,
    trlist,
    gfdir,
    return_gfs=False,
    residuals=False,
    OnlyGetFullGF=False,
    max_t_d=200,
):
    '''
    Perform the actual W-phase inversion.

    :param float t_d: Time delay
    :param tuple cmtloc: (lat, lon, depth) of the centroid's location.
    :param tuple periods: (Ta,Tb), passband periods.
    :param array MRF: Moment rate function.
    :param array observed_displacements: Array containing concatenated traces of observed disp.
    :param array trlen: Array containing the length of each trace.
    :param dict metadata: Dictionary with the metadata of each station.
    :param list trlist: List with the station id which will contribute to the inv.
    :param gfdir: Path to the greens functions or a GreensFunctions instances.
    :param bool hdf5_flag: *True* if the greens functions are stored in HDF5, *False* if they are in SAC.
    :param bool return_gfs: *True* if the greens functions should be returned.
    :param bool residuals: *True*, return the 'raw' misfit from the least squares inversion,
        *False*, return the 'relative' misfit as a percentage: i.e. the 'raw' misfit divided by
        the norm of the synthetics.
    :param bool OnlyGetFullGF: *True*, return the greens function matrix for the maximum time delay (*max_t_d*)
        without performing the inversion, *False*, perform the inversion for the time delay given
        by *t_d*.
    :param numeric max_t_d: Maximum time delay to consider if *OnlyGetFullGF = True*.

    :return: What is returned depends on the values of the parameters as described below.

        - If *OnlyGetFullGF = True*, then just the greens function matrix.
        - Otherwise a tuple containing:

            0. moment tensor components in Nm ['RR', 'PP', 'TT', 'TP', 'RT', 'RP']
            #. misfit (percent L2 norm misfit error of the solution), and

            If *return_gfs = True*

            2. The greens function matrix also.
    '''
    logger.debug("core_inversion(t_d=%d, c=(%.2f, %.2f, %.2f), #tr=%d)", t_d, cmtloc[0], cmtloc[1], cmtloc[2], len(trlist))
    if isinstance(gfdir, GreensFunctions):
        greens = gfdir
    else:
        greens = GreensFunctions(gfdir)
    delta = greens.delta

    hyplat, hyplon, hypdep = cmtloc
    Ta, Tb = periods

    # Index of the first value that will be valid after convolution with MRF.
    first_valid_time = int(len(MRF)/2.)

    # the indices of the beginning and end of each trace in observed displacements
    indexes =  np.array(np.concatenate((np.array([0.]), np.cumsum(trlen))), dtype='int')


    if OnlyGetFullGF:
        max_t_d = int(max_t_d)
        Nst = len(trlen)
        GFmatrix = np.zeros((np.array(trlen, dtype=int).sum() + max_t_d*Nst, 5))
        tb = 0
    else:
        GFmatrix = np.zeros((np.array(trlen, dtype=int).sum(), 5))

    #### Inversion:
    for i, trid in enumerate(trlist):
        trmeta = metadata[trid]
        trlat = trmeta['latitude']
        trlon = trmeta['longitude']

        try:
            sta  = trmeta['sta']
        except KeyError:
            sta = None

        dist = locations2degrees(hyplat, hyplon, trlat, trlon)
        distm, azi, bazi  =  gps2dist_azimuth(hyplat, hyplon, trlat, trlon)
        t_p =  trmeta.get('ptime')

        if not t_p:
            from obspy.taup.taup import getTravelTimes
            t_p =  getTravelTimes(dist,hypdep)[0]['time']

        # Select greens function and perform moment tensor rotation in the azimuth
        # DETAIL: The greens functions are stored for zero azimuth (which we can do
        # because of "symetry properties of a sphere". Then to recover the correct
        # moment tensor (or equivalently synthetic trace) we rotate them. For details
        # see Kanamori and Rivera 2008.
        azi_r = -azi*np.pi/180
        synth = greens.select_rotated(trid[-1], dist, hypdep, azi_r)

        ##################################
        # Extract W phase from synthetics.

        # Use only what we think is noise to calculate and compensate the mean
        synth -= np.mean(synth[..., :60], axis=-1)[:, None]

        # Convolve with moment rate function
        synth = ndimage.convolve1d(synth, MRF, axis=-1, mode='nearest')\
            [:, first_valid_time-2:-first_valid_time]

        # Pad data back to original length
        fillvec1 = synth[:, 0, None] * np.ones(first_valid_time)
        fillvec2 = synth[:, -1, None] * np.ones(first_valid_time)
        synth = np.concatenate((fillvec1, synth, fillvec2), axis=-1)
        synth -= np.mean(synth[:, :60], axis=-1)[:, None]

        # Bandpass filter to extract W phase
        def wphase_filter(trace):
            return bandpassfilter(trace, delta, 4, 1./Tb, 1./Ta)
        synth = np.apply_along_axis(wphase_filter, -1, synth)

        if OnlyGetFullGF:
            synth = ltrim(synth, t_p - max_t_d, delta)
            synth = synth[:, :trlen[i] + max_t_d]
        else:
            synth = ltrim(synth, t_p - t_d, delta)
            synth = synth[:, :trlen[i]]

        if OnlyGetFullGF:
            ta = tb
            tb = ta + max_t_d + (indexes[i+1] - indexes[i])
        else:
            ta = indexes[i]
            tb = indexes[i+1]

        # the first of the following lines are due to constraint that volume does not change
        GFmatrix[ta:tb, 0] = synth[0][:] - synth[1][:]
        GFmatrix[ta:tb, 1] = synth[2][:] - synth[1][:]
        GFmatrix[ta:tb, 2] = synth[4][:]
        GFmatrix[ta:tb, 3] = synth[5][:]
        GFmatrix[ta:tb, 4] = synth[3][:]

    if OnlyGetFullGF:
        return GFmatrix

    # perform the inversion
    inversion = lstsq(GFmatrix, observed_displacements, rcond=None)
    M = inversion[0]

    # construct the synthetics
    syn = GFmatrix.dot(M)

    # set the 'sixth' element of the moment tensor (corresponds to the
    # constraint added 20 lines or so above)
    M = np.insert(M, 2, -M[0]-M[1])

    # extract the residuals and scale if required
    if residuals:
        misfit = float(inversion[1])
    else:
        misfit = 100.*np.sqrt(np.sum((syn-observed_displacements)**2)/np.sum(syn**2))

    if return_gfs:
        return M, misfit, GFmatrix

    return M, misfit



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



def core_inversion_wrapper(allarg):
    '''
    Wrapper for :py:func:`core_inversion` for use with
    :py:class:`multiprocessing.Pool`.

    :return: None
    '''

    if isinstance(allarg[-1], dict):
        kwargs = allarg[-1]
        arg = allarg[:-1]
    else:
        kwargs = {}
        arg = allarg

    try:
        return core_inversion(*arg,**kwargs)
    except Exception:
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))


def minimize_misfit(inputs, processes=None):
    """Given a list of input tuples for core_inversion, find the one with the
    lowest misfit.

    :param inputs: an iterable of tuples to provided to core_inversion as args.
    :param processes: number of parallel processes to use.
    :return: A tuple (index_of_minimizer, moment_tensor, misfit, results list)."""
    with ProcessPoolExecutor(max_workers=processes) as pool:
        results = list(pool.map(core_inversion_wrapper, inputs))
    misfits = np.array([x[1] for x in results])
    i_min = misfits.argmin()
    return i_min, results[i_min][0], results[i_min][1], results


def remove_individual_traces(tol, M, GFmatrix, observed_displacements,
                                  trlist, trlen):
    '''
    Remove any trace with an indivudual misfit higher than *tol*.

    :param float tol: The maximum acceptable value for a misfit.
    :param array M: Moment tensor.
    :param array GFmatrix: Greens function matrix.
    :param array Observed Disp: The observed displacements.
    :param list trlist: ???
    :param list trlen: List containing the lengths of the traces in *trlist*.

    :return: Tuple containing the

        0. greens function matrix,
        #. observed displacements,
        #. trace list, and
        #. list of the lengths of the traces in the trace list

        for the traces which remain.
    '''

    GFmatrix_fix = np.zeros((1,5))
    ObservedDisp_fix = np.zeros(1)
    trlist_fix = []
    trlen_fix = []
    M = np.delete(M,2)

    # construct the synthetic trace
    syn = GFmatrix.dot(M)
    j = 0
    for i, trid in enumerate(trlist):
        # get the observed and synthetic traces for a channel
        n1 = j
        n2 = trlen[i] + j
        obs = observed_displacements[n1:n2]
        syn2 = syn[n1:n2]

        # calculate the misfit for the channel
        misfit = 100.*np.linalg.norm(syn2-obs)\
                 /np.linalg.norm(syn2)

        j += trlen[i]

        # add the trace to the remaining traces if misfit is within tolerance.
        if misfit <= tol :
            ObservedDisp_fix = np.append(ObservedDisp_fix, obs)
            GFmatrix_fix = np.append(GFmatrix_fix, GFmatrix[n1:n2,:], axis=0)
            trlist_fix.append(trid)
            trlen_fix.append(trlen[i])

    # build the result
    GFmatrix = GFmatrix_fix[1:,:]
    observed_displacements = ObservedDisp_fix[1:]
    trlist = trlist_fix[:]
    trlen = trlen_fix[:]

    return GFmatrix, observed_displacements, trlist, trlen

def get_timedelay_misfit_wrapper(args):
    return get_timedelay_misfit(*args)



def get_timedelay_misfit(t_d, GFmatrix, trlen, observed_displacements, max_t_d):
    try:
        max_t_d = int(max_t_d)
        t_d = int(t_d)
        t_d2 = max_t_d - t_d
        GFmatrix_sm = np.zeros((np.array(trlen,dtype=int).sum(),5))
        cumtrlens = np.concatenate(([0], np.array(trlen,dtype=int).cumsum()))
        tb_fm = 0
        for i_ntr, ta in enumerate(cumtrlens[:-1]):
            l = trlen[i_ntr]
            tb = ta + l
            ta_fm = tb_fm
            tb_fm = ta_fm +  max_t_d + tb-ta
            GFmatrix_sm[ta:tb] = GFmatrix[ta_fm + t_d2 :ta_fm + t_d2 + l]

        inversion = lstsq(GFmatrix_sm, observed_displacements, rcond=None)
        return inversion[1][0]
    except Exception:
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))


def make_result(cls, M: Sequence[float], misfit: float, **extra):
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
    logger.info("%s Result:" % cls.__name__)
    logger.info("%5s % 10s % 10s % 10s", "", "r", "t", "p")
    logger.info("%5s %+.3e %+.3e %+.3e", "r", M[0], M[3], M[4])
    logger.info("%5s % 10s %+.3e %+.3e", "t", "",   M[1], M[5])
    logger.info("%5s % 10s % 10s %+.3e", "p", "",   "",   M[2])
    logger.info("misfit: %.0f%%", misfit)
    M2 = np.asarray(M) ** 2
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
