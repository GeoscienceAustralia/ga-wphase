# -*- coding: utf-8 -*-
import os
import logging
import numpy as np
import pandas as pd
from typing import Optional
from collections import defaultdict
from traceback import format_exc

# to avoid: Exception _tkinter.TclError
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
try:
    from obspy.imaging.beachball import (beachball as plot_beachball,
                                         aux_plane, mt2plane, MomentTensor)
except ImportError:
    from obspy.imaging.beachball import (Beachball as plot_beachball,
                                         AuxPlane as aux_plane,
                                         MT2Plane as mt2plane,
                                         MomentTensor)


from wphase.psi.seismoutils import get_azimuths, azimuthal_gap
from wphase.plotting import (
    plot_grid_search,
    plot_preliminary_fit,
    plot_screening_stage,
    plot_station_coverage,
    plot_waveforms,
)
from wphase.psi import model
from wphase import settings

logger = logging.getLogger(__name__)


try:
    # If we can import pyinstrument imported, then profile
    from pyinstrument import Profiler

    class WPInvProfiler(object):
        """WPInvProfiler."""

        def __init__(self, working_dir=None):
            """__init__.

            Parameters
            ----------
            working_dir :
                working_dir
            """
            self.html = None
            self.working_dir = working_dir

        def __enter__(self):
            """__enter__."""
            self.profiler = Profiler()  # or Profiler(use_signal=False), see below
            self.profiler.start()

        def __exit__(self, exc_type, esc_value, traceback):
            """__exit__.

            Parameters
            ----------
            exc_type :
                exc_type
            esc_value :
                esc_value
            traceback :
                traceback
            """
            self.profiler.stop()
            self.html = self.profiler.output_html()
            if self.working_dir is not None:
                with open(
                    os.path.join(self.working_dir, "timings.html"), "w"
                ) as timings_file:
                    timings_file.write(self.profiler.output_html())


except Exception:
    import cProfile, pstats, io

    class WPInvProfiler(object):
        """WPInvProfiler."""

        def __init__(self, *args, **kwargs):
            """__init__.

            Parameters
            ----------
            args :
                args
            kwargs :
                kwargs
            """
            self.html = None
            self.sort_by = "cumulative"  #'tottime'

        def __enter__(self):
            """__enter__."""
            self.profiler = cProfile.Profile()
            self.profiler.enable()

        def __exit__(self, exc_type, esc_value, traceback):
            """__exit__.

            Parameters
            ----------
            exc_type :
                exc_type
            esc_value :
                esc_value
            traceback :
                traceback
            """
            self.profiler.disable()
            s = io.StringIO()
            ps = pstats.Stats(self.profiler, stream=s).sort_stats(self.sort_by)
            ps.print_stats()
            self.html = "<pre>{}</pre>".format(s.getvalue())


class NoProfiler(object):
    def __init__(self, *args, **kwargs):
        self.html = None

    def __enter__(self):
        pass

    def __exit__(self, exc_type, esc_value, traceback):
        pass


class OutputDict(defaultdict):
    """
    A dict with an additional method for adding warnings and which converts any
    numpy array in itself or the items it holds (recursively) to a python list.
    """

    def __init__(self):
        super(OutputDict, self).__init__(OutputDict)

    def __setitem__(self, key, value):
        def _recurse(d):
            if isinstance(d, dict) and not isinstance(d, OutputDict):
                res = OutputDict()
                for k, v in d.items():
                    super(OutputDict, res).__setitem__(k, _recurse(v))
                return res
            elif isinstance(d, np.ndarray):
                return list(d)
            return d

        return super(OutputDict, self).__setitem__(key, _recurse(value))

    def add_warning(self, msg):
        """
        Add a the warning message given by *msg*.

        Warnings are accumulated in a list under the key "Warnings"
        """

        if "Warnings" not in self:
            self["Warnings"] = []
        self["Warnings"].append(msg)

    def as_dict(self, item=None):
        if item is None:
            return {k: None if v is None else self.as_dict(v) for k, v in self.items()}
        if isinstance(item, OutputDict):
            return item.as_dict()
        return item


def convert_to_antelope(M, cenloc):
    """
    Convert Roberto's wphase output to antelope's moment tensor format.
    Plus calculate a few values based on the moment tensor.

    :param M: Moment tensor
    :param cenloc: The centroid location, (cenlat,cenlon,cendep)
    :return: antelope's moment tensor info in a dictionary.
    """
    M2 = np.asarray(M)**2
    m0 = np.sqrt(0.5*(M2[0]+M2[1]+M2[2])+M2[3]+M2[4]+M2[5])
    mag = 2./3.*(np.log10(m0)-9.10)

    np1 = mt2plane(MomentTensor(M, 0))
    np2 = aux_plane(np1.strike, np1.dip, np1.rake)

    results = model.AntelopeMomentTensor(
        tmpp=M[2],
        tmrp=M[4],
        tmrr=M[0],
        tmrt=M[3],
        tmtp=M[5],
        tmtt=M[1],

        scm=m0,
        drmag=mag,
        drmagt='Mww',

        drlat=cenloc[0],
        drlon=cenloc[1],
        drdepth=cenloc[2],

        str1=np1.strike,
        dip1=np1.dip,
        rake1=np1.rake,

        str2=np2[0],
        dip2=np2[1],
        rake2=np2[2],

        auth=settings.AUTHORITY,
    )

    try:
        results.dc, results.clvd = decomposeMT(M)
    except Exception:
        import traceback
        logger.warning("Error computing DC/CLVD decomposition: %s",
                       traceback.format_exc())

    return results


def plot_and_save_beachball(M, working_dir, OL):
    beachBallPrefix = os.path.join(
        working_dir,
        settings.BEACHBALL_PREFIX
    )
    plot_beachball(
        M,
        width=400,
        outfile=beachBallPrefix + "_OL{}.png".format(OL),
        format='png'
    )
    plt.close('all') # obspy doesn't clean up after itself...


def post_process_wpinv(
    output: model.WPhaseResult,
    WPOL,
    working_dir,
    eqinfo: model.Event,
    metadata,
    make_maps=True,
    make_plots=True):
    traces = None
    MT = None
    mtResult: Optional[model.OL2Result] = output.OL3 or output.OL2

    if output.OL1:
        traces = output.OL1.used_traces
        prelim = output.OL1.preliminary_calc_details
        fname = os.path.join(working_dir, settings.PRELIM_FIT_PREFIX) + ".png"
        plot_preliminary_fit(eqinfo, filename=fname, **prelim)
    else:
        logger.warning("Could not find preliminary calculation details in result.")

    if mtResult:
        MT = mtResult.moment_tensor
        traces = mtResult.used_traces
        output.QualityParams = model.Quality(
            azimuthal_gap=azimuthal_gap(
                get_azimuths(metadata, traces, (eqinfo.latitude, eqinfo.longitude))
            ),
            number_of_stations=len(set(trid.split(".")[1] for trid in traces)),
            number_of_channels=len(traces),
        )

        if make_plots and len(traces) > 0:
            output.WphaseResultPlots = plot_waveforms(
                working_dir,
                settings.RESULTS_TRACES_PREFIX,
                mtResult.synthetic_displacements,
                mtResult.observed_displacements,
                mtResult.trace_lengths,
            )

        else:
            output.add_warning(
                "Could not create wphase waveforms plot - no traces to plot!"
            )

    if make_plots and output.OL2 is not None:
        try:
            plot_and_save_beachball(output.OL2.moment_tensor, working_dir, OL=2)
        except Exception:
            output.add_warning(
                "Failed to create beachball for OL2. {}".format(format_exc())
            )

    if make_plots and output.OL3 is not None:
        try:
            plot_and_save_beachball(output.OL3.moment_tensor, working_dir, OL=3)
        except Exception:
            output.add_warning(
                "Failed to create beachball for OL3. {}".format(format_exc())
            )

    if isinstance(mtResult, model.OL3Result):
        cenloc = mtResult.centroid

        output.MomentTensor = convert_to_antelope(MT, cenloc)

        # Only 3 has cenloc...
        output.Centroid = model.CentroidLocation(
            depth=round(cenloc[2], 1),
            latitude=round(cenloc[0], 3),
            longitude=round(cenloc[1], 3),
            time=mtResult.centroid_time,
        )

        if make_maps:
            # draw the grid search plot
            coords = np.asarray(mtResult.grid_search_candidates)
            misfits = np.array([x[1] for x in mtResult.grid_search_results])
            N_grid = len(misfits)
            lats, lons, depths = coords.T
            depths_unique  = sorted(set(depths))
            N_depths = len(depths_unique)
            misfits_depth_mat = np.zeros((int(N_grid/N_depths),N_depths))
            latlon_depth_mat = np.zeros((int(N_grid/N_depths),2,N_depths))
            ##We will sum the misfits over the depths
            for i_col,depth in enumerate(depths_unique):
                i_depth = np.where(depths == depth)
                misfits_depth_mat[:,i_col] = misfits[i_depth]
                latlon_depth_mat[:,:,i_col] = coords[i_depth,:2]

            ##This should be the same for all depths
            latlon_depth_grid =  latlon_depth_mat[:,:,0]
            #Suming all the depths
            misfits_depth_mat =  misfits_depth_mat.sum(axis=1)
            scaled_field = misfits_depth_mat/misfits_depth_mat.min()

            gridSearchPrefix = os.path.join(working_dir, settings.GRID_SEARCH_PREFIX)
            plot_grid_search(
                (eqinfo.longitude, eqinfo.latitude),
                (cenloc[1], cenloc[0]),
                latlon_depth_grid,
                scaled_field,
                s=100./scaled_field**2,
                c=scaled_field,
                zorder=999,
                filename=gridSearchPrefix
            )

    if traces and make_maps:
        try:
            # Make a plot of the station distribution
            lats = [metadata[trid]['latitude'] for trid in traces]
            lons = [metadata[trid]['longitude'] for trid in traces]
            stationDistPrefix = os.path.join(
                working_dir,
                settings.STATION_DISTRIBUTION_PREFIX)
            plot_station_coverage(
                (eqinfo.latitude, eqinfo.longitude),
                lats,
                lons,
                mt=MT,
                filename=stationDistPrefix + '.png')
        except Exception:
            output.add_warning("Failed to create station distribution plot. {}".format(format_exc()))

    if make_maps:
        for i, stage in enumerate(output.ScreeningStages):
            outfile = f"{working_dir}/stage{i:02d}_{stage.name}.png"
            logger.info(f"Creating screening map {outfile}...")
            def station_row(pair):
                trid, passed = pair
                try:
                    return {
                        "lat": metadata[trid]['latitude'],
                        "lon": metadata[trid]['longitude'],
                        "station": trid.rsplit(".", 2)[0],
                        "passed": passed,
                    }
                except Exception:
                    return None
            rows = map(station_row, stage.station_results.items())
            plot_screening_stage(
                (eqinfo.latitude, eqinfo.longitude),
                pd.DataFrame.from_records([row for row in rows if row]),
                mt=stage.mt,
                filename=outfile,
                title=stage.info or stage.name,
            )


def decomposeMT(M):
    """Given a deviatoric (i.e. trace-free) moment tensor (specified as a
    6-element list in the CMT convention as usual), compute the percentage
    double-couple and compensated linear vector dipole components.

    Written following this paper:
    Vavryčuk, V. Moment tensor decompositions revisited. J Seismol 19, 231–252
    (2015). https://doi.org/10.1007/s10950-014-9463-y

    :returns: A tuple ``(DC, CLVD)`` of relative scale factors between 0 and 1.
    """
    mt = MomentTensor(M, 0)
    eigs, _ = np.linalg.eig(mt.mt)
    M1, M2, M3 = np.sort(eigs)[::-1] # M1 >= M2 >= M3

    # Since we're working with deviatoric moment tensors, we are assuming
    # M_ISO=0 and we don't have to worry about the destinction between the
    # Silver&Jordan and Knopoff&Randall decompositions.
    M_CLVD = (2./3.)*(M1 + M3 - 2*M2)
    M_DC = (1./2.)*(M1 - M3 - abs(M1 + M3 - 2*M2))
    M = abs(M_CLVD) + abs(M_DC)
    return abs(M_DC)/M, abs(M_CLVD)/M
