from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import range
from builtins import object
import os
import logging
import numpy as np
from collections import defaultdict

# to avoid: Exception _tkinter.TclError
import matplotlib
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
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


from wphase.psi import seismoutils
from wphase.psi.plotutils import plot_field, stacov, make_preliminary_fit_plot
from wphase.psi.model import OL1, OL2, OL3
from wphase import settings

logger = logging.getLogger(__name__)


if settings.PROFILE_WPHASE:
    try:
        # If we can import pyinstrument imported, then profile
        from pyinstrument import Profiler
        class WPInvProfiler(object):
            def __init__(self, wphase_output, working_dir):
                self.wphase_output = wphase_output
                self.working_dir = working_dir

            def __enter__(self):
                self.profiler = Profiler() # or Profiler(use_signal=False), see below
                self.profiler.start()

            def __exit__(self, exc_type, esc_value, traceback):
                self.profiler.stop()
                self.wphase_output[settings.WPINV_PROFILE_OUTPUT_KEY] = \
                    self.profiler.output_html()
                if self.working_dir is not None:
                    with open(os.path.join(self.working_dir, 'timings.html'), 'w') as timings_file:
                        timings_file.write(self.profiler.output_html())

    except Exception:
        import cProfile, pstats, io
        class WPInvProfiler(object):
            def __init__(self, wphase_output, *args, **kwargs):
                self.wphase_output = wphase_output
                self.sort_by = 'cumulative'#'tottime'

            def __enter__(self):
                self.profiler = cProfile.Profile()
                self.profiler.enable()

            def __exit__(self, exc_type, esc_value, traceback):
                self.profiler.disable()
                s = io.StringIO()
                ps = pstats.Stats(self.profiler, stream=s).sort_stats(self.sort_by)
                ps.print_stats()
                self.wphase_output[settings.WPINV_PROFILE_OUTPUT_KEY] = {
                    'css':'',
                    'js':'',
                    'body':'<pre>{}</pre>'.format(s.getvalue())}
else:
    class WPInvProfiler(object):
        def __init__(self, *args, **kwargs): pass
        def __enter__(self): pass
        def __exit__(self, exc_type, esc_value, traceback): pass





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

        Warnings are accumulated in a list under the key
        :py:data:`wphase.settings.WPHASE_WARNINGS_KEY`.
        """

        if settings.WPHASE_WARNINGS_KEY not in self:
            self[settings.WPHASE_WARNINGS_KEY] = []
        self[settings.WPHASE_WARNINGS_KEY].append(msg)

    def as_dict(self, item=None):
        if item is None:
            return {k: None if v is None else self.as_dict(v) for k, v in self.items()}
        if isinstance(item, OutputDict):
            return item.as_dict()
        return item





def wpinv_for_eatws(M, cenloc):
    """
    Convert Roberto's wphase output to antelope's moment tensor format.
    Plus calculate a few values based on the moment tensor.

    :param M: Moment tensor
    :param cenloc: The centroid location, (cenlat,cenlon,cendep)
    :return: antelope's moment tensor info in a dictionary.
    """

    results = {}
    results['tmpp'] = M[2]
    results['tmrp'] = M[4]
    results['tmrr'] = M[0]
    results['tmrt'] = M[3]
    results['tmtp'] = M[5]
    results['tmtt'] = M[1]

    # from roberto's code
    M2 = M*M
    m0 = np.sqrt(0.5*(M2[0]+M2[1]+M2[2])+M2[3]+M2[4]+M2[5])
    mag = 2./3.*(np.log10(m0)-9.10)

    results['scm'] = m0
    results['drmag'] = mag
    results['drmagt'] = 'Mww'

    results['drlat'] = cenloc[0]
    results['drlon'] = cenloc[1]
    results['drdepth'] = cenloc[2]

    moment_tensor = MomentTensor(M, 0)
    nodalplane = mt2plane(moment_tensor)
    results['str1'] = nodalplane.strike
    results['dip1'] = nodalplane.dip
    results['rake1'] = nodalplane.rake

    np2 = aux_plane(
        nodalplane.strike,
        nodalplane.dip,
        nodalplane.rake)
    results['str2'] = np2[0]
    results['dip2'] = np2[1]
    results['rake2'] = np2[2]

    results['auth'] = settings.GA_AUTHORITY

    return results





def post_process_wpinv(
    res,
    wphase_output,
    WPOL,
    working_dir,
    eqinfo,
    metadata,
    make_maps=True,
    make_plots=True):

    prelim = res.preliminary_calc_details

    if prelim:
        fname = os.path.join(working_dir, settings.WPHASE_PRELIM_FIT_PREFIX) + '.png'
        make_preliminary_fit_plot(eqinfo, filename=fname, **prelim)
    else:
        logger.warning("Could not find preliminary calculation details in result.")

    M_OL2 = None
    WPOL = 1
    # extract the results to local variables
    if isinstance(res, OL2):
        WPOL = 2
        M = res.moment_tensor
        obs = res.observed_displacements
        syn = res.synthetic_displacements
        traces = res.trace_lengths
    else:
        traces = list(res.preliminary_calc_details['trids'])

    if isinstance(res, OL3):
        WPOL = 3
        M_OL2 = wphase_output['OL2'].pop('M')
        cenloc = res.centroid
        inputs_latlon = res.grid_search_candidates
        moments = res.grid_search_results

        if make_plots:
            try:
                # Display the beachball for OL2
                beachBallPrefix = os.path.join(
                    working_dir,
                    settings.WPHASE_BEACHBALL_PREFIX)
                plot_beachball(M_OL2, width=400,
                    outfile = beachBallPrefix + "_OL2.png", format='png')
                plt.close('all') # obspy doesn't clean up after itself...
            except Exception:
                wphase_output.add_warning("Failed to create beachball for OL2.")

    if 'OL2' in wphase_output:
        wphase_output['OL2'].pop('M', None)

    if WPOL >= 2:
        wphase_output['QualityParams']['azimuthal_gap'] = seismoutils.AzimuthalGap(
                metadata,
                traces,
                (eqinfo['lat'], eqinfo['lon']))[0]
        wphase_output['QualityParams']['number_of_stations'] = len(set(
            trid.split('.')[1] for trid in traces))
        wphase_output['QualityParams']['number_of_channels'] = len(traces)

    if make_plots:
        try:
            # Display the beachball for the output level achieved
            beachBallPrefix = os.path.join(working_dir, "{}_OL{}".format(
                settings.WPHASE_BEACHBALL_PREFIX, WPOL))
            plot_beachball(M, width=400,
                outfile = beachBallPrefix + ".png", format='png')
            plt.close('all') # obspy doesn't clean up after itself...
        except Exception:
            wphase_output.add_warning("Failed to create beachball for OL{}.".format(
                WPOL))

    if make_maps:
        try:
            # Make a plot of the station distribution
            hyplat =  eqinfo['lat']
            hyplon =  eqinfo['lon']
            lats = [metadata[trid]['latitude'] for trid in traces]
            lons = [metadata[trid]['longitude'] for trid in traces]
            stationDistPrefix = os.path.join(
                working_dir,
                settings.WPHASE_STATION_DISTRIBUTION_PREFIX)
            stacov(
                (hyplat,hyplon),
                lats,
                lons,
                mt=M,
                filename=stationDistPrefix + '.png')
        except Exception:
            wphase_output.add_warning("Failed to create station distribtuion plot.")

    if WPOL >= 2 and make_plots and len(traces):
        # Secondly the wphase traces plot, syn Vs obs
        class PlotContext(object):
            """
            Handle a single plot.
            """
            def __init__(self, name):
                self.name = name
                self.xticks = []
                self.xlabel = []

            def add_tick(self, tick_pos, label):
                self.xticks.append(tick_pos)
                self.xlabel.append(label)

            def save_image(self, folder, syn, obs, formats=['png']):
                fig = matplotlib.figure.Figure(figsize=(14, 7))
                FigureCanvas(fig)
                ax = fig.add_subplot(1, 1, 1)
                ax.set_title("Wphase Results (Red: Synthetic, Blue: Observed)")
                ax.plot(syn, color="red")
                ax.plot(obs, color="blue")
                ax.set_xticks(self.xticks)
                ax.set_xticklabels(self.xlabel, rotation=90, fontsize='xx-small')
                ax.set_xlim((0, len(obs)))
                for offset in self.xticks:
                    ax.axvline(offset, color='0.1')
                for fmt in formats:
                    fig.savefig(
                        os.path.join(folder, '{}.{}'.format(self.name, fmt)),
                        dpi=120,
                        format=fmt,
                        bbox_inches='tight')

        class CreatePlots(object):
            """
            Generate the wphase result plots.

            This is achieved by constructing an instance of this class.

            :param folder: The folder to write the plots to.
            :param prefix: Prefix for file names.
            :param syn: The synthetic data (for all stations).
            :param obs: The observed data (for all stations).
            :param n_traces: The number of traces (i.e. stations).

            .. note:: This uses *traces* which is
                currently defined in the enclosing scope.
            """
            def __init__(self, folder, prefix, syn, obs, n_traces):
                self.folder = folder
                self.prefix = prefix
                self.syn = syn
                self.obs = obs
                self.n_traces = n_traces
                self.all_traces = PlotContext(prefix)
                self.images = []
                self.images.append(((1, n_traces), prefix + '.png'))
                self.n_traces_in_curr = 0
                if n_traces > settings.N_TRACES_PER_RESULT_PLOT:
                    self.n_subplots_done = 0
                    self.start_next_plot(0)
                else:
                    self.curr_sub_plot = None

                # generate the plots
                offset = 0
                for trid, trlen in traces.items():
                    offset += trlen
                    self.add_tick(offset, trid.split('.')[1])

            def add_tick(self, tick_pos, label):
                """
                Add a tick (which corresponds to a station).
                """
                self.all_traces.add_tick(tick_pos, label)
                if self.curr_sub_plot is not None:
                    self.n_traces_in_curr += 1
                    self.curr_sub_plot.add_tick(tick_pos - self.start_index, label)
                    if self.n_traces_in_curr == settings.N_TRACES_PER_RESULT_PLOT:
                        self.save_curr_subplot(tick_pos)
                        self.start_next_plot(tick_pos)

            def save_curr_subplot(self, end_index):
                """
                Save the plot for a subset of traces.
                """
                slc = slice(self.start_index, end_index)
                self.curr_sub_plot.save_image(self.folder, self.syn[slc], self.obs[slc])
                self.images.append((self.next_plot_range, self.curr_sub_plot.name + '.png'))
                self.n_subplots_done += 1

            def start_next_plot(self, end_index):
                """
                Start a plot for the next subset of traces.
                """
                self.n_traces_in_curr = 0
                self.start_index = end_index
                first = self.n_subplots_done * settings.N_TRACES_PER_RESULT_PLOT + 1
                last = (self.n_subplots_done + 1) * settings.N_TRACES_PER_RESULT_PLOT
                last = min(last, self.n_traces)
                self.next_plot_range = (first, last)#'{} to {}'.format(first, last)
                self.curr_sub_plot = PlotContext('{}_{}_{}'.format(self.prefix, first, last))

            def __del__(self):
                """
                Save the plots of all traces and the last set of traces.
                """
                self.all_traces.save_image(self.folder, self.syn, self.obs, ['png'])
                if self.curr_sub_plot is not None and self.n_traces_in_curr:
                    self.save_curr_subplot(len(self.syn))
                wphase_output[settings.RESULTS_PLOTS_KEY] = self.images

        CreatePlots(
            working_dir,
            settings.WPHASE_RESULTS_TRACES_PREFIX,
            syn,
            obs,
            len(traces))

    elif make_plots:
        wphase_output.add_warning('Could not create wphase results plot.')

    if WPOL==3:
        results = wpinv_for_eatws(M, cenloc)
        wphase_output['MomentTensor'] = results

        # Only 3 has cenloc...
        wphase_output['Centroid'] = {}
        wphase_output['Centroid']['depth'] = round(cenloc[2],1)
        wphase_output['Centroid']['latitude'] = round(cenloc[0],3)
        wphase_output['Centroid']['longitude'] = round(cenloc[1],3)

        if make_maps:
            # draw the grid search plot
            inputs = inputs_latlon
            N_grid = len(inputs)
            misfits = np.array([moments[i_grid][1] for i_grid in range(N_grid)])
            coords = np.array([inputs[i_grid][2] for i_grid in range(N_grid)])
            lats, lons, depths = coords[:,:].T
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

            gridSearchPrefix = os.path.join(working_dir, settings.WPHASE_GRID_SEARCH_PREFIX)
            plot_field((eqinfo['lon'], eqinfo['lat']),
                       (cenloc[1], cenloc[0]),
                       latlon_depth_grid,
                       scaled_field,
                       s=100./scaled_field**2,
                       c=scaled_field,
                       zorder=999,
                       filename=gridSearchPrefix)
    else:
        results = None

    return results
