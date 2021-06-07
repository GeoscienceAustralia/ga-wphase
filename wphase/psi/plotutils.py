# -*- coding: utf-8 -*-
'''
Some utilities to make nice plots.
'''
from __future__ import absolute_import, print_function

from builtins import range

import logging
from collections import OrderedDict

import numpy as np
import matplotlib.colors as mcols
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from scipy.interpolate import griddata
from obspy.imaging.beachball import beach

try:
    # cartopy is an optional dependency
    from cartopy import crs
    from cartopy.feature import NaturalEarthFeature
except ImportError:
    pass

logger = logging.getLogger(__name__)

def get_boundaries(lats,lons, xpnd=0.7):
    lnmn = lons.min()
    lnmx = lons.max()
    ltmn = lats.min()
    ltmx = lats.max()
    lllon = 0.5*(lnmn+lnmx)-xpnd*(lnmx-lnmn)
    lllat = 0.5*(ltmn+ltmx)-xpnd*(ltmx-ltmn)
    urlon = 0.5*(lnmn+lnmx)+xpnd*(lnmx-lnmn)
    urlat = 0.5*(ltmn+ltmx)+xpnd*(ltmx-ltmn)
    return lllon, urlon, lllat, urlat



def plot_field(
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
    Plot the locations considered in the grid search.

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
        fig = Figure()
        FigureCanvas(fig)

    lats, lons = latlons.T
    boundaries = get_boundaries(lats, lons, xpnd = 0.7)

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

def stacov(location, lats, lons, mt=None, filename=None, fig=None):
    """Plot a station coverage map."""

    (elat, elon) = location
    if not fig:
        fig = Figure()
        FigureCanvas(fig)

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

def make_preliminary_fit_plot(eqinfo,
                              strike, average_amplitude, anisotropy,
                              corrected_amplitudes, azimuths,
                              filename, **kwargs):
    fig = Figure()
    FigureCanvas(fig)
    ax = fig.add_subplot(111)
    title = 'Preliminary Magnitude Fit'
    try:
        title += ' for %s' % eqinfo['id']
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
