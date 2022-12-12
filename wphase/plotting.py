# -*- coding: utf-8 -*-
"""Some utilities to make nice plots."""
from __future__ import absolute_import, print_function
import os

from builtins import range

import logging
from collections import OrderedDict
from cartopy.crs import PlateCarree

import numpy as np
import matplotlib.colors as mcols
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from scipy.interpolate import griddata
from obspy.imaging.beachball import beach

from wphase import settings

try:
    # cartopy is an optional dependency
    from cartopy import crs
    from cartopy.feature import NaturalEarthFeature
except ImportError:
    pass

logger = logging.getLogger(__name__)


def make_figure(*args, **kwargs):
    """Create a new matplotlib figure with an in-memory aggregate canvas."""
    fig = Figure(*args, **kwargs)
    FigureCanvasAgg(fig)
    return fig

def get_boundaries(lats, lons, padding=1.4):
    """Given a set of points (specified as the two iterables ``lats`` and
    ``lons``), return a bounding box containing all of them."""
    lnmn = lons.min()
    lnmx = lons.max()
    ltmn = lats.min()
    ltmx = lats.max()
    lllon = 0.5*(lnmn+lnmx)-0.5*padding*(lnmx-lnmn)
    lllat = 0.5*(ltmn+ltmx)-0.5*padding*(ltmx-ltmn)
    urlon = 0.5*(lnmn+lnmx)+0.5*padding*(lnmx-lnmn)
    urlat = 0.5*(ltmn+ltmx)+0.5*padding*(ltmx-ltmn)
    return lllon, urlon, lllat, urlat



def plot_grid_search(
    location,
    centroid,
    latlons,
    field,
    filename,
    plot_type='scatter',
    show_lats=True,
    show_lons=True,
    fig=None,
    clevs=None,
    **kwargs):
    '''
    Plot a map showing the locations considered in the grid search.

    :param location: (lon, lat) coords of the earthquake.
    :param latlons: Tuples contain the lat/lon pairs.
    :type latlons: Two column numpy array.
    :param field: The field to plot. Must have the dimensions commensurate to *latlons*.
    :type field: Two dimensional numpy array.
    :param str plot_type: Either "scatter" (for points) or "contour" (for contours).
    :param bool show_lats: Should latitudes be shown on the y-axis.
    :param bool show_lons: Should longitudes be shown on the x-axis.
    :param ax: Axes for figure to add the plot to.
    :type ax: None or :py:class:`matplotlib.Axes`.
    :param list clevs: The levels for the contour plot. Ignored if *plot_type* is "scatter".
    :param str topfile: Path to file to load 'etopo' topography from.
    '''

    if not fig:
        fig = make_figure()

    lats, lons = latlons.T
    boundaries = get_boundaries(lats, lons)

    proj = crs.Mercator(central_longitude=centroid[0])
    coords = crs.PlateCarree()

    ax = fig.add_axes((0.05, 0.18, 0.95, 0.75), projection=proj)
    ax.set_extent(boundaries)

    ocean = NaturalEarthFeature('physical', 'ocean', '110m')
    land = NaturalEarthFeature('physical', 'land', '110m')
    countries = NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '110m')

    ax.add_feature(ocean, facecolor='aqua')
    ax.add_feature(land, facecolor='#FF9900')
    ax.add_feature(countries, edgecolor='grey', facecolor='none')
    ax.coastlines(resolution='110m')
    ax.gridlines(crs=coords, draw_labels=True,
                 linestyle=(0, [5,5]), color='grey')

    #Plotting field
    if plot_type == 'scatter':
        field_plot = ax.scatter(lons, lats, transform=coords, **kwargs)
    else:
        # contour maps are no longer used, so I haven't bothered reimplementing
        # them in cartopy yet.
        raise NotImplementedError()
    fig.colorbar(field_plot, orientation='vertical')
    ax.scatter(*location, c='y', s=1000, marker="*", zorder=1000, transform=coords)
    ax.scatter(*centroid, c='w', s=1000, marker="*", zorder=1001, transform=coords)
    grid_legend="\n".join((
        "Colorbar indicates normalized centroid misfit (1 is minimum)",
        "Yellow star: Hypocenter location",
        "White star: optimal centroid location",
    ))
    fig.text(.5, -.05, grid_legend, horizontalalignment='center')
    fig.savefig(filename, bbox_inches='tight')

def plot_station_coverage(location, trids, lats, lons, mt=None, filename=None, fig=None):
    """Plot a map showing the stations used in a W-Phase solution."""

    (elat, elon) = location
    if not fig:
        fig = make_figure()

    proj = crs.Orthographic(central_longitude=elon, central_latitude=elat)
    coords = crs.PlateCarree()
    ax = fig.add_subplot(111, projection=proj)
    ax.set_global()

    ax.coastlines(resolution='110m')
    land = NaturalEarthFeature('physical', 'land', '110m')
    ocean = NaturalEarthFeature('physical', 'ocean', '110m')
    ax.add_feature(land, facecolor='#FF9900')
    ax.add_feature(ocean, facecolor='aqua')
    ax.gridlines(crs=coords, draw_labels=False,
                 color='grey', linestyle=(0, [5,5]), linewidth=0.4)

    ax.scatter(lons, lats, transform=coords,
               s=60, c='blue', marker='o', edgecolors='none', zorder=10)

    latlon_crs = PlateCarree()
    for lon, lat, trid in zip(lons, lats, trids):
        xy = ax.projection.transform_point(lon, lat, src_crs=latlon_crs)
        ax.annotate(trid, xy)

    from warnings import warn
    if mt is not None:
        projected_xy = proj.transform_point(elon, elat, coords)
        bball = beach(mt, xy=(0,0), linewidth=1, facecolor='r', zorder=10000,
                      width=8e5)
        ax.add_collection(bball)
    else:
        logger.warning("No MT provided; station coverage map will have no beachball")
        ax.plot(elon, elat, 'y*', markersize=40, transform=coords)

    if filename:
        fig.savefig(filename, dpi=100, transparent=True, bbox_inches='tight')
        logger.info("Wrote %s successfully", filename)

def plot_preliminary_fit(eqinfo, strike, average_amplitude, anisotropy,
                         corrected_amplitudes, azimuths, filename, **kwargs):
    """Plot amplitude against azimuth to visualize the preliminary magnitude
    fit."""
    fig = make_figure()
    ax = fig.add_subplot(111)
    title = 'Preliminary Magnitude Fit'
    try:
        title += ' for %s' % eqinfo.id
    except Exception:
        pass
    ax.set_title(title)
    ax.set_xlabel(u'Azimuth (°)')
    ax.set_ylabel('P2P Amplitude (mm, corrected for attenuation)')
    ax.set_xlim([0, 360])

    ax.scatter(azimuths, 1e3*np.asarray(corrected_amplitudes), label='station amplitudes',
               c='red', marker='^', alpha=0.3)

    rule_style = dict(linestyle='--', alpha=0.2, lw=1)
    hlabel = u'fitted P2P ampl = %.1emm' % (2e3*average_amplitude)
    vlabel = u'fitted strike = %.0f°' % strike
    ax.axhline(2e3*average_amplitude, label=hlabel, **rule_style)
    ax.axvline(strike % 360, label=vlabel, **rule_style)
    ax.legend()

    x = np.arange(0, 360)
    y = 1e3*(2*average_amplitude - anisotropy*np.cos(2*np.pi*(x-strike)/180))
    ax.plot(x, y, color='blue', label='best fit')

    fig.savefig(filename, dpi=100, bbox_inches='tight')
    logger.info("Wrote %s successfully", filename)

class WaveformPlot(object):
    def __init__(self, name):
        """Helper class for visualizing synthetic vs observed waveforms."""
        self.name = name
        self.xticks = []
        self.xlabel = []

    def add_tick(self, tick_pos, label):
        self.xticks.append(tick_pos)
        self.xlabel.append(label)

    def save_image(self, folder, syn, obs, formats=['png']):
        fig = make_figure(figsize=(14, 7))
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


class MultiWaveformPlotter(object):
    def __init__(self, folder, prefix, syn, obs, traces):
        """
        Helper class used by plot_waveforms.

        :param folder: The folder to write the plots to.
        :param prefix: Prefix for file names.
        :param syn: The synthetic data (for all stations).
        :param obs: The observed data (for all stations).
        :param traces:
            An OrderedDict mapping station names to the number of
            samples in the corresponding waveforms.
        """
        self.folder = folder
        self.prefix = prefix
        self.syn = syn
        self.obs = obs
        self.traces = traces
        self.n_traces = len(traces)
        self.all_traces = WaveformPlot(prefix)
        self.images = []
        self.n_traces_in_curr = 0
        self.n_subplots_done = 0
        self.curr_sub_plot = None

    def plot(self):
        self.images.append(((1, self.n_traces), self.prefix + ".png"))
        if self.n_traces > settings.N_TRACES_PER_RESULT_PLOT:
            self.start_next_plot(0)

        # generate the plots
        offset = 0
        for trid, trlen in self.traces.items():
            offset += trlen
            self.add_tick(offset, trid.split('.')[1])

        self.all_traces.save_image(self.folder, self.syn, self.obs, ["png"])
        if self.curr_sub_plot is not None and self.n_traces_in_curr:
            self.save_curr_subplot(len(self.syn))

        return self.images

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
        if self.curr_sub_plot:
            slc = slice(self.start_index, end_index)
            self.curr_sub_plot.save_image(self.folder, self.syn[slc], self.obs[slc])
            self.images.append((self.next_plot_range, self.curr_sub_plot.name + ".png"))
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
        self.curr_sub_plot = WaveformPlot('{}_{}_{}'.format(self.prefix, first, last))


def plot_waveforms(folder, prefix, syn, obs, traces):
    """Plot synthetic vs observed waveforms.

    Generates one plot including all waveforms, but also a sequence of plots
    each showing the waveforms for N_TRACES_PER_RESULT_PLOT stations.

    :param folder: The folder to write the plots to.
    :param prefix: Prefix for file names.
    :param syn: The synthetic data (for all stations).
    :param obs: The observed data (for all stations).
    :param traces:
        An OrderedDict mapping station names to the number of
        samples in the corresponding waveforms."""
    return MultiWaveformPlotter(folder, prefix, syn, obs, traces).plot()
