''' This module contains useful fuctions to handle seismic data
    with obspy
'''
from __future__ import absolute_import, print_function

from future import standard_library
standard_library.install_aliases()
from builtins import range
import urllib.request, urllib.error, urllib.parse, sys
import logging
import math
import numpy as np

logger = logging.getLogger(__name__)

try:
    from obspy.geodetics import locations2degrees
except ImportError:
    from obspy.core.util.geodetics import locations2degrees

try:
    from obspy.geodetics.base import gps2dist_azimuth
except ImportError:
    from obspy.core.util.geodetics import gps2DistAzimuth as gps2dist_azimuth



def get_azimuths(meta, trlist, loc):
    '''
    Return the azimuths in degrees given metadata dict, list of
    station ids and location of the epicentre.
    '''
    return [gps2dist_azimuth(loc[0], loc[1], meta[t]['latitude'], meta[t]['longitude'])
            for t in trlist]

def station_pruning(meta, trlist, cutoffs=[1.,2.,5], units='deg'):
    '''
    Given a list of id trlist and its respective metadata dictionary meta, it
    will return a reduced list containing only stations at a certain minimum distance
    from each other. The cutoff list, with units given by this keyword arg, will
    be iteratively checked to remove a station from a pair closer than cutoff[i].
    If unit is km the result should be more precise (but slower) owing to the use
    of gps2dist_azimuth from obspy.

    Note: trlist should not cointain duplicated ids.
    '''

    red_trlist = list(trlist) #copy of trlist

    for cutoff in cutoffs:
        for i_trid, trid1 in enumerate(trlist[:-1]):
            for j_trid,trid2 in enumerate(trlist[i_trid + 1:]):
                lat1, lon1  = meta[trid1]['latitude'], meta[trid1]['longitude']
                lat2, lon2  = meta[trid2]['latitude'], meta[trid2]['longitude']
                if units == 'deg':
                    dist = locations2degrees(lat1,lon1,lat2,lon2)
                elif units == 'km':
                    dist = 1.e-3*gps2dist_azimuth(lat1,lon1,lat2,lon2)[0]


                if dist <= cutoff:
                    try:
                        red_trlist.remove(trid1)
                    except Exception:
                        pass
        trlist = list(red_trlist)# copying reduced list
    return  trlist


def azimuthal_gap(azis):
    """Compute the (largest) azimuthal gap of a list of azimuths.

    Parameters
    ----------
    azis : Iterable[float]
        Station azimuths in degrees.
    """
    if len(azis) < 2:
        return 360.
    azis = np.asarray(azis) % 360
    azis.sort()
    gaps = (azis - np.roll(azis, 1)) % 360
    return gaps.max()
