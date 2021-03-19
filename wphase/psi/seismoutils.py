''' This module contains useful fuctions to handle seismic data
    with obspy
'''
from __future__ import print_function
from __future__ import absolute_import

from future import standard_library
standard_library.install_aliases()
from builtins import range
import urllib.request, urllib.error, urllib.parse, sys
import math as mat
import numpy as np

from wphase import logger

try:
    from obspy.geodetics import locations2degrees
except ImportError:
    from obspy.core.util.geodetics import locations2degrees

try:
    from obspy.geodetics.base import gps2dist_azimuth
except ImportError:
    from obspy.core.util.geodetics import gps2DistAzimuth as gps2dist_azimuth

from . import datautils as DU

def getMetadataFromSeed(sp, trid, datetime=None):
    '''
    Given a parser object from a dataless seed and a station id , this function returns a
    dictionary containing every station id and: Latitude, Longitude, Type of Transfer function,
    Elevation, Azimuth, Dip, poles, zeros, gain (normalization factor) and sensitivity.

    PS: This function is based on getPAZ and getCoordinates of OBSPY but it adds more
    information to the output. Note that transfer function is specified.
    '''

    blks = sp._select(trid,datetime)
    poles = None
    for blk in blks:
        if blk.id == 58:
                if blk.stage_sequence_number == 0:
                    sens = blk.sensitivity_gain
                elif blk.stage_sequence_number == 1:
                    seisgain = blk.sensitivity_gain
                elif blk.stage_sequence_number == 2:
                    diggain = blk.sensitivity_gain



        if blk.id == 52:
            lat = blk.latitude
            lon = blk.longitude
            Elev= blk.elevation
            azi = blk.azimuth
            dip = blk.dip

        if blk.id == 53 and poles == None:
            tf = blk.transfer_function_types
            gain = blk.A0_normalization_factor
            poles = []
            for i in range(blk.number_of_complex_poles):
                try:
                    p = complex(blk.real_pole[i], blk.imaginary_pole[i])
                except TypeError:
                    p = complex(blk.real_pole, blk.imaginary_pole)
                poles.append(p)

            zeros = []
            for i in range(blk.number_of_complex_zeros):
                try:
                    z = complex(blk.real_zero[i], blk.imaginary_zero[i])
                except TypeError:
                    z = complex(blk.real_zero, blk.imaginary_zero)
                zeros.append(z)




    metadata   =  { 'latitude'       : lat,
                     'longitude'      : lon,
                     'elevation'      : Elev,
                     'transfer_function': tf,
                     'azimuth'        : azi,
                     'dip'            : dip,
                     'zeros'          : zeros,
                     'poles'          : poles,
                     'gain'           : gain,
                     'sensitivity'    : sens
 #                   'seismometer_gain': seisgain,
  #                   'digitizer_gain' : diggain
                   }

    return metadata


def hypo2dist(hyploc1,hyploc2,r=6371.):
    '''
    Given 2 tuples with hypcenter coordinates (lat, lon, depth) returns the
    distance in KM between them.
    '''

    r1 = r-hyploc1[2]
    r2 = r-hyploc2[2]
    torad = np.pi/180.
    #dellon = np.abs(torad*(hyploc1[1]-hyploc2[1]))
    theta = torad*locations2degrees(hyploc1[0],hyploc1[1],\
                                         hyploc2[0],hyploc2[1])
    #theta2 = np.arccos(np.sin(torad*hyploc1[0])*np.sin(torad*hyploc2[0])
    #               + np.cos(torad*hyploc1[0])*np.cos(torad*hyploc2[0])*np.cos(dellon))
    ##Law of Cosines:
    l = np.sqrt(r1**2.+r2**2.-2.*r2*r1*np.cos(theta))
    #l2 = np.sqrt(r1**2.+r2**2.-2.*r2*r1*np.cos(theta2))
    return l#,l2



def AzimuthalGap(META,trlist, location):
    '''
    Return the maximun Azimuthal Gap given metadata dict, list of
    station ids and location of the epicentre.
    '''
    (eplat,eplon) = location
    if len(trlist) < 2:
        return 359.999, []
    azis = np.empty(len(trlist) + 1)
    for i, trid in enumerate(trlist):
        trMETA = META[trid]
        trlat = trMETA['latitude']
        trlon = trMETA['longitude']
        azi =  gps2dist_azimuth(eplat, eplon, trlat, trlon)[1]
        azis[i] = azi
    azis[-1] = azis[:-1].min() +360.0
    args = np.argsort(azis)
    trlist = [trlist[i] for i in args[:-1]]
    trlist.append(trlist[0])
    azis.sort()
    azis_diff = azis[1:] - azis[:-1]
    gap =  azis_diff.max()
    azi_cov = []
    for i in range(len(trlist[:-1])):
        #key = trlist[i] + ' -> ' + trlist[i+1]
        #azi_cov[key] = azis_diff[i]
        azi_cov.append([trlist[i], trlist[i+1], azis_diff[i]])
    return  gap, azi_cov


def GetVirtualNetworkStations(VN):
    '''
    Given a Virtual Network (e.g. "FDSN") it returns
    a list with the stations in the Network.
    '''

    url = "http://service.iris.edu/irisws/virtualnetwork/1/query?code=_" + VN +"&format=CSV"
    #url = "http://www.iris.edu/vnets?vnet=_" + VN +"&vout=CSV"
    lines = urllib.request.urlopen(url).read().split("\n")
    stations = []
    for line in lines:
        if line.startswith('_'):
            fields = line.rstrip().split(',')
            net = fields[1]
            sta = fields[2]
            stations.append(net + "." + sta)

    return stations

def resample_Ntraces(Ntraces, DecFac):
    '''
    This function recomputes the proper length of each trace (given as a
    Ntraces list) after we applied a decimation to the full array.
    '''

    Ntraces_dec = np.empty(len(Ntraces))
    res = 0.
    new_val = Ntraces[0]
    if np.any(np.array(Ntraces) <= DecFac):
        logger.warning("At least one trace is smaller that the Decimation Factor")
    for i in range(len(Ntraces)):
        if not i == 0:
            res = round((mat.ceil(new_val/DecFac)
                        - new_val/DecFac)*DecFac)
        new_val = Ntraces[i]-res
        Ntraces_dec[i] = mat.ceil(new_val/DecFac)
    return Ntraces_dec
    ##

def rot_12_NE(st, META):
    '''
    Performs a  12 -> NE rotation. META is a dictionary with the metadata
    and st is a  stream.
    '''

    st2 = st.select(channel="??1")
    for tr in st2:
        id1 = tr.id
        id2 = id1[:-1] + '2'
        tr1 = tr
        try:
            tr2 = st.select(id=id2)[0]
        except IndexError:
            st.remove(tr)
            logger.warning(tr.id, "Channel 2 not found. Impossible to rotate")
            continue
        timeA = max(tr1.stats.starttime,tr2.stats.starttime)
        timeB = min(tr1.stats.endtime,tr2.stats.endtime)
        tr1.trim(timeA,timeB)
        tr2.trim(timeA,timeB)
        azi = META[id1]['azimuth']
        tr1.data, tr2.data = DU.Rot2D(tr1.data,tr2.data,-azi)
    return st

def station_pruning(META,trlist, cutoffs=[1.,2.,5], units='deg'):
    '''
    Given a list of id trlist and its respective metadata dictionary META, it
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
                lat1, lon1  = META[trid1]['latitude'], META[trid1]['longitude']
                lat2, lon2  = META[trid2]['latitude'], META[trid2]['longitude']
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
    return   trlist


def get_new_coordinates(lat1, lon1, azi, dist):
    '''computes the location (lat,lon) of a new point given the a location and a
       distance, azimuth pair
       UPDATE: I think this is not a right approach. The Azimimuth will not project
               as I thought. Better to use  geographiclib.
    '''

    r = 6371.
    dist_rad = np.pi/180.*dist
    azi_rad = np.pi/180.*azi
    lat1_rad = np.pi/180.*lat1
    lon1_rad = np.pi/180.*lon1

    dp = np.sqrt(2*r**2*(1-np.cos(dist_rad)))
    y = dp*np.cos(azi_rad)
    x = dp*np.sin(azi_rad)

    dlat = np.arccos(1-(y**2)/(2*r**2))
    dlon = np.arccos(1-(x**2)/(2*r**2))
    if azi > 90 and azi < 270:
        dlat = -dlat
    if azi > 180 and azi < 360:
        dlon = -dlon

    lat2_rad = dlat + lat1_rad
    lon2_rad = dlon + lon1_rad

    return 180/np.pi*lat2_rad, 180/np.pi*lon2_rad

