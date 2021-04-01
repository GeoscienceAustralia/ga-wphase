"""Core implementation of the W-Phase inversion."""
from __future__ import absolute_import, print_function

from builtins import str
from builtins import range
import sys, os, glob, traceback
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Pool
from timeit import default_timer as timer
import numpy as np
from numpy.linalg import lstsq
from scipy.interpolate import interp1d
from scipy.signal import triang
from scipy import ndimage
from scipy.signal import detrend
import h5py
from obspy.core import read, Stream

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

from wphase import logger
from .greens import GreensFunctions
from .bandpass import bandpassfilter

class WPInvWarning(Exception):
    """
    A 'terminal' condition that we do not consider to be an error.

    This should be raised if there is a problem with in input data or some other
    issue from which we cannot recover.

    :param str msg: The reason for termination.
    """

    def __init__(self, msg):
        super(WPInvWarning, self).__init__(msg)

class RTdeconvError(Exception):
    pass

def wpinv(
    st,
    metadata,
    eqINFO,
    gfdir,
    OL=1,
    processes=None,
    output_dic=None):
    '''
    This function is the main function that will compute the inversion. For
    details of the the algorithm and logic, see `here
    <https://pubs.er.usgs.gov/publication/70045140>`_.

    :param st: containing the raw data.
    :type st: :py:class:`obspy.core.stream.Stream`

    :param dict metadata: Station metadata. Each key is a station ID and
        the values are a dictionaries with the metadata. Each dictionary
        should look like\:

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
                sum(Ntrace) = len(ObservedDisp) = len(syn).

        - OL = 3, Same as OL2 plus:
            5. Optimal centroid location (lat, lon, dep).
            #. Time delay/Half duration of the source (secs).
            #. latitude-longitude grid used to find the optimal centroid
                location
            #. Inversion result for each grid point in latlons
            #. Dictionary with relavant information about the data processing
                so it can redone easily outside this function.
    '''

    if output_dic is None:
        output_dic = {}

    # utility for adding warnings to output_dict
    def add_warning(msg):
        if type(output_dic) is dict:
            # note that the key must be the same as atws_wphase.settings.WARNINGS_KEY,
            # but we can't import from that, here, and don't want to import
            # from any non-standard package in atws_wphase.settings, so we just
            # have to manually keep these in sync.
            WARNINGS_KEY = 'Warnings'
            if WARNINGS_KEY not in output_dic:
                output_dic[WARNINGS_KEY] = []
            output_dic[WARNINGS_KEY].append(msg)
        else:
            output_dic.add_warning(msg)

    # used to reject a channel if the fitting of the instrument response is
    # greater than this.  see "Source inversion of W phase - speeding up
    # seismic tsunami warning - Kanamori - 2008" in the doc folder of
    # atws_wphase for details on the fitting of the instrument response.
    response_misfit_tol = 5. #Percent

    # time cutoff, seconds/degree from t_p for the Wphase window.
    wphase_cutoff = 15

    # traces with peak-to-peak amplitude outside of the range
    # [median_rejection_coeff[0] * m, median_rejection_coeff[1] * m], where m
    # is the median peak-to-peak amplitude over all stations, will be rejected
    median_rejection_coeff = [0.1, 3]

    # Minimum number of stations required to trigger the inversion.
    N_st_min = 4

    # Maximum time delay in seconds. We'll search for the optimal t_d up to
    # this value.
    max_t_d = 200.

    # stations where the misfits (100*sqrt(sum(synthetic-observed)^2 /
    # sum(synthetic)^2)) are greater than misfit_tol[0] will be rejected. If
    # any stations remain, then those with misfits greater than misfit_tol[1]
    # will be rejected, and so on.
    misfit_tol = [300, 200, 100]#, 80, 60, 50]# 40, 30]

    # Number of points in which the finer search will be performed in OL3. That
    # is the Nfinersearch lat/lons with the lowest misfits will have a finer
    # scale search done in their neighbourhoods. This comes from earlier
    # versions of the code and may not work at present.
    Nfinersearch = 0

    # After removing bad traces need this # of channels to continue
    minium_num_channels = 10

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

    hyplat = eqINFO['lat']
    hyplon = eqINFO['lon']
    hypdep = eqINFO['dep']
    orig = eqINFO['time']

    # Deal with extremly shallow preliminary hypocenter
    if hypdep < 10.:
        hypdep = 10.

    #In case of multiple locations we favor
    # '00'. This may be improved if we select the one with the longer period
    # sensor.
    st_sel = st.select(location = '00')
    st_sel+= st.select(location = '--')
    st_sel+= st.select(location = '') #Check this.
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
        t2 = t1 + dist*wphase_cutoff

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
        (om0, h, G) = getCOEFFfit(trmeta['sensitivity'],freq,AmpfromPAZ)
        if not np.all(np.isfinite([om0, h, G])):
            logger.warning("Imposible to get Coeff. Skipping %s", tr.id)
            trlist.remove(tr.id)
            continue

        AmpfromCOEFF= np.abs(omega*omega / \
                (omega*omega + 2j*h*om0*omega - om0*om0))

        # L2 norm:
        misfit = 100*np.linalg.norm(AmpfromPAZ-AmpfromCOEFF) \
                / np.linalg.norm(AmpfromPAZ)

        if misfit >  response_misfit_tol:
            logger.warning('Bad fitting for response. Skipping:\n{}\t {: E}\n'.format(
                    tr.id,
                    misfit))
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
    median_p2p = np.median(tr_p2p)
    accepted_traces = [i for i in range(len(tr_p2p))
                         if tr_p2p[i] < median_rejection_coeff[1]*median_p2p
                         and tr_p2p[i] > median_rejection_coeff[0]*median_p2p]
    if not EnoughAzimuthalCoverage():
        raise WPInvWarning("Lack of azimuthal coverage. Aborting.")
    if len(accepted_traces) < N_st_min:
        raise WPInvWarning("Lack of stations. Aborting.")

    logger.info("Traces accepted for preliminary magnitude calculation: {}"
                .format(len(accepted_traces)))

    # Get the preliminary mag:
    tr_p2p_con = [tr_p2p[i] for i in accepted_traces]
    DIST_con   = [DIST[i] for i in accepted_traces]
    AZI_con    = [AZI[i] for i in accepted_traces]
    trlist_pre_con = [trlist_pre[i] for i in accepted_traces]
    pre_param = preliminarymagnitude(tr_p2p_con, DIST_con, AZI_con)

    pre_wp_mag = pre_param[0]
    pre_M0 =  pre_param[1]
    pre_strike = pre_param[2]
    t_h = pre_param[3]

    if pre_wp_mag < 6.5:
        add_warning("Preliminary magnitude less than 6.5... setting it to 6.5 and continuing")
        pre_wp_mag = 6.5

    logger.info("OL1:")
    logger.info("Preliminary W-phase magnitude for the event: %.7f", pre_wp_mag)

    output_dic['OL1'] = {}
    output_dic['OL1']['magnitude'] = round(pre_wp_mag,1)
    output_dic['OL1']['nstations'] = len(accepted_traces)
    output_dic['OL1']['used_traces'] = trlist

    if OL==1:
        return pre_wp_mag, trlist_pre_con, tr_p2p_con

    #############  Output Level 2    #######################################
    #############  Moment Tensor based on preliminary hypercenter (PDE) ####
    # Much of what follows is the same as what was done above, but we will be
    # using a different set of stations and can handle multiple components per
    # station (i.e. horizontals also)
    ########################################################################

    # Redefine and define some values according to the pre_wp_mag
    Ta, Tb =   GetCornerFreqsFromMag(pre_wp_mag)
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
        t2 = t1 + dist*wphase_cutoff
        tr = st_sel.select(id = trid)[0]

        tr.data = np.array(tr.data,dtype=float)
        if trmeta['transfer_function'] == "B":
            AmpfromPAZ  = Vpaz2freq(trmeta,freq/2./np.pi)  # hz to rad*hz
        elif trmeta['transfer_function'] == "A":
            AmpfromPAZ  = Vpaz2freq(trmeta,freq)
        else:
            logger.warning("Unknown transfer function. Skipping %s", tr.id)
            trlist.remove(tr.id)
            continue

        (om0, h, G) = getCOEFFfit(trmeta['sensitivity'],freq,AmpfromPAZ)
        if not np.all(np.isfinite([om0, h, G])):
            logger.warning("Impossible to get Coeff. Skipping %s", tr.id)
            trlist.remove(tr.id)
            continue

        DATA_INFO[tr.id] = [om0, h, G, dt, t1, t2]

        try:
            tr.data, coeff = RTdeconv(
                tr,
                om0,
                h,
                G,
                dt,
                corners=4,
                baselinelen=60./dt,
                taperlen= 10.,
                fmin = 1./Tb,
                fmax = 1./Ta,
                get_coef = True)

        except RTdeconvError as e:
            logger.warning("Error deconvolving trace %s: %s", tr.id, str(e))
            trlist.remove(tr.id)
            continue

        #Check length of the trace:
        tr.trim(t1,t2)
        if len(tr)== 0:
            "Empty trace. Skipping", tr.id
            trlist.remove(tr.id)
            continue
        tr_p2p.append(tr[:].max()-tr[:].min())

        DIST.append(dist)

    DIST = np.array(DIST)
    trlist = [trlist[i] for i in np.argsort(DIST)]
    tr_p2p = [tr_p2p[i] for i in np.argsort(DIST)]

    # Rejecting outliers:
    ObservedDisp = np.array([]) # observed displacements vector.
    Ntrace = [] # A list with the length of each station data.
    trlist2 = trlist[:]

    median_AMP = np.median(tr_p2p)
    for i, amp in enumerate(tr_p2p):
        if (amp > median_AMP*median_rejection_coeff[1]
        or amp < median_AMP*median_rejection_coeff[0]):
            trlist2.remove(trlist[i])

        else:
            tr = st_sel.select(id = trlist[i])[0]
            ObservedDisp =  np.append(ObservedDisp, tr.data[:],0)
            Ntrace.append(len(tr))
    trlist = trlist2[:]

    'number of traces for OL2 {}'.format(len(trlist))

    #### Inversion:
    # Search for the optimal t_d
    ##I need to test the optimal way to do this. map is taking too much time

    # time delays for which we will run inversions (finally choosing the one with
    # lowest misfit)
    time_delays = np.arange(1., max_t_d)

    # extract the greens function matrix for all time delays. Note that this does not
    # perform an inversion, because OnlyGetFullGF=True.
    GFmatrix = core_inversion(
        t_h,
        0,
        (hyplat, hyplon, hypdep),
        orig,
        (Ta, Tb),
        MRF,
        ObservedDisp,
        Ntrace,
        metadata,
        trlist,
        gfdir=gfdir,
        OnlyGetFullGF=True,
        Max_t_d=max_t_d,
    )

    # inputs for the TIME DELAY search
    inputs = [(t_d_test, GFmatrix, Ntrace, ObservedDisp, max_t_d) for t_d_test in time_delays]
    with ProcessPoolExecutor() as pool:
        misfits = list(pool.map(get_timedelay_misfit_wrapper,inputs))

    # Set t_d (time delay) and and t_h (half duration) to optimal values:
    mis_min =  np.array(misfits).argmin()
    output_dic['misfits'] = {'min':mis_min, 'array':misfits}
    t_d = t_h = time_delays[mis_min]
    MRF = MomentRateFunction(t_h, dt)
    logger.info("Source time function, time delay: %d, %f", len(MRF), t_d)

    #### Removing individual bad fitting. this recursively removes stations with misfits
    # outside of the acceptable range defined by the variable misfit_tol, which was
    # described above.
    M, misfit, GFmatrix = core_inversion(
        t_h,
        t_d,
        (hyplat, hyplon, hypdep),
        orig,
        (Ta, Tb),
        MRF,
        ObservedDisp,
        Ntrace,
        metadata,
        trlist,
        gfdir=gfdir,
        GFs=True)

    # Remove bad traces
    for tol in misfit_tol:
        # Make sure there are enough channels
        GFmatrix, ObservedDisp, trlist, Ntrace = remove_individual_traces(
            tol,
            M,
            GFmatrix,
            ObservedDisp,
            trlist,
            Ntrace)

        if len(trlist) < minium_num_channels:
            output_dic.pop('OL1', None)
            msg = "Only {} channels with possibly acceptable fits. Aborting.".format(len(trlist))
            logger.error(msg)
            raise WPInvWarning(msg)

        M, misfit, GFmatrix = core_inversion(
            t_h,
            t_d,
            (hyplat, hyplon, hypdep),
            orig,
            (Ta,Tb),
            MRF,
            ObservedDisp,
            Ntrace,
            metadata,
            trlist,
            gfdir=gfdir,
            GFs=True
        )

    syn = (M[0]*GFmatrix[:,0] + M[1]*GFmatrix[:,1] +
           M[3]*GFmatrix[:,2] + M[4]*GFmatrix[:,3] +
           M[5]*GFmatrix[:,4])

    logger.info("OL2:")

    output_dic['OL2'] = MT_result(M, misfit, hypdep, t_d)

    if OL==2:
        return  M, ObservedDisp, syn, trlist, Ntrace, GFmatrix, DATA_INFO
    else:
        output_dic['OL2']['M'] = M

    if len(Ntrace) == 0:
        add_warning("Could not calculate OL3: no data within tolerance")
        return  M, ObservedDisp, syn, trlist, Ntrace, GFmatrix, DATA_INFO

    #############  Output Level 3   #############################
    ###Moment Tensor based on grid search  ######################
    ###for the optimal centroid's location #####################

    logger.info("building latlon grid")
    lat_grid, lon_grid = get_latlon_for_grid(hyplat, hyplon, dist_lat=3.0,
                                             dist_lon=3.0, delta=0.8)
    logger.debug("Grid size: %d * %d", len(lon_grid), len(lat_grid))
    latlons = [(lat, lon) for lat in lat_grid for lon in lon_grid]

    inputs = [(t_h, t_d, (lat, lon, hypdep), orig, (Ta, Tb), MRF,
              ObservedDisp, Ntrace, metadata, trlist, greens,
              True, dict(residuals=False))
              for lat, lon in latlons]

    with ProcessPoolExecutor() as pool:
        latlon_search = list(pool.map(core_inversion_wrapper, inputs))

    misfits_latlon = np.array([latlon_search[i][1] for i in range(len(latlons))])
    cenlat, cenlon = latlons[misfits_latlon.argmin()]
    moments = latlon_search  #Compatibility

    deps_grid = get_depths_for_grid(hypdep, greens)
    logger.debug("Depth grid size: %d", len(deps_grid))

    inputs = [(t_h, t_d, (cenlat, cenlon, depth), orig, (Ta,Tb),
               MRF, ObservedDisp, Ntrace, metadata, trlist, greens,
               True, dict(residuals=False))
              for depth in deps_grid]

    with ProcessPoolExecutor() as pool:
        depth_search = list(pool.map(core_inversion_wrapper, inputs))

    misfits_depth = np.array([depth_search[i][1] for i in range(len(deps_grid))])
    cendep = deps_grid[misfits_depth.argmin()]


    ###Final inversion!!

    M, misfit, GFmatrix = core_inversion(t_h, t_d, (cenlat, cenlon, cendep),
                                         orig, (Ta,Tb), MRF, ObservedDisp,
                                         Ntrace, metadata, trlist, gfdir=gfdir,
                                         GFs=True)

    syn = (M[0]*GFmatrix[:,0] + M[1]*GFmatrix[:,1]
          + M[3]*GFmatrix[:,2] + M[4]*GFmatrix[:,3]
          + M[5]*GFmatrix[:,4])


    logger.info("OL3:")
    output_dic['OL3'] = MT_result(M, misfit, cendep, t_d)

    cenloc = (cenlat,cenlon,cendep)
    return M, ObservedDisp, syn, trlist, Ntrace, cenloc, t_d, latlons, moments, DATA_INFO


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


def GetCornerFreqsFromMag(Mw):
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



def preliminarymagnitude(tr_p2p, dists, azis):
    '''
    Compute the preliminary magnitude.

    :param list tr_p2p: Peak to peak amplitudes for each trace.
    :param list dists: station - epi distances in deg for each trace.
    :param list azis: station - epi azimuths in degrees for each trace.

    :return: tuple containing:

        0. Preliminary W-phase Mw magnitude.
        #. Scalar moment in dyne * cm.
        #. Prelimary estimation for the strike in degrees.
        #. Preliminary half duration for the moment rate func in seconds.
    '''

    # Parameters:
    # Attenuation relation for W-phase:
    distance = np.array([5., 10., 20., 30., 40., 50.,
                         60., 70., 80., 90.5])
    atten = np.array([1.5, 1.4, 1.2, 1.0, 0.7, 0.56,
                       0.61, 0.56, 0.50 , 0.52 ])

    # Empirical coeff. for W-phase mag. log A = m*(M_w) + n
    m = np.log10(0.9/0.065)
    n = np.log10(0.9)-m*8.5

    azis = np.array(azis)
    N = len(tr_p2p)
    azis *= np.pi/180.
    atten_func =  interp1d(distance,atten, kind = 'cubic')

    ## We will solve the system Mx = B in the least squares sense.
    B = tr_p2p/atten_func(dists)
    M = np.zeros((N,3))
    M[:,0] = 1
    M[:,1] = np.cos(2*azis)
    M[:,2] = np.sin(2*azis)
    x = lstsq(M,B, rcond=None)[0]
    amp = x[0]/2.

    ###We need to add 3 because we use meters instead of mm.
    pre_wp_mag = (np.log10(amp)-n+3)/m
    #Scalar moment in dyne * cm
    M0 = 10 ** (1.5 * pre_wp_mag + 16.1)
    t_h = 1.2  * 10** -8. * M0 **(1./3.)

    # I need to check  this. I do not know the convention
    # for the strike I should use.
    pre_strike = 0.5*np.arctan2(x[2],x[1])*180./np.pi

    return  pre_wp_mag, M0, pre_strike, t_h



def EnoughAzimuthalCoverage():
    ##### TO BE IMPLEMENTED
    return True



def getCOEFFfit(sens, freq, resp):
    """
    Fit the amplitude of the instrument response at a given set of frequencies.

    :param float sens: The instrument's sensitivity.
    :param array freq: The target frequencies.
    :param resp: The amplitudes of the response in freq.
    :type resp: array of the same length as *freq*.

    :return: Tuple containing the three time domain deconvolution coefficients.
    """

    X = (2.*np.pi*freq)**2
    Y = X*X*(1./resp/resp-1.)
    fit = np.polyfit(X,Y,1)
    G = sens
    om0 = fit[1]**(.25)
    h = np.sqrt(fit[0]/4./np.sqrt(fit[1]) + 0.5)
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


def core_inversion(t_h, t_d, cmtloc,orig, periods, MRF,
                   ObservedDisp, Ntrace, metadata, trlist,
                   gfdir,
                   GFs = False, residuals = False,
                   OnlyGetFullGF=False, Max_t_d=200,
                   ):
    '''
    Perform the actual W-phase inversion.

    :param float t_h: Half duration of the MRF in seconds.
    :param float t_d: Time delay
    :param tuple cmtloc: (lat, lon, depth) of the centroid's location.
    :param orig: Origin time.
    :type orig: :py:class:`obspy.core.obspy.core.utcdatetime.UTCDateTime`
    :param tuple periods: (Ta,Tb), passband periods.
    :param array MRF: Moment rate function.
    :param array ObservedDisp: Array containing concatenated traces of observed disp.
    :param array Ntraces: Array containing the length of each trace.
    :param dict metadata: Dictionary with the metadata of each station.
    :param list trlist: List with the station id which will contribute to the inv.
    :param gfdir: Path to the greens functions or a GreensFunctions instances.
    :param bool hdf5_flag: *True* if the greens functions are stored in HDF5, *False* if they are in SAC.
    :param bool GFs: *True* if the greens functions should be returned.
    :param bool residuals: *True*, return the 'raw' misfit from the least squares inversion,
        *False*, return the 'relative' misfit as a percentage: i.e. the 'raw' misfit divided by
        the norm of the synthetics.
    :param bool OnlyGetFullGF: *True*, return the greens function matrix for the maximum time delay (*Max_t_d*)
        without performing the inversion, *False*, perform the inversion for the time delay given
        by *t_d*.
    :param numeric Max_t_d: Maximum time delay to consider if *OnlyGetFullGF = True*.

    :return: What is returned depends on the values of the parameters as described below.

        - If *OnlyGetFullGF = True*, then just the greens function matrix.
        - Otherwise a tuple containing:

            0. moment tensor components in Nm ['RR', 'PP', 'TT', 'TP', 'RT', 'RP']
            #. misfit (percent L2 norm misfit error of the solution), and

            If *GFs = True*

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

    # the indicies of the beginning and end of each trace in observed displacements (ObservedDisp)
    indexes =  np.array(np.concatenate((np.array([0.]), np.cumsum(Ntrace))), dtype='int')


    if OnlyGetFullGF:
        Max_t_d = int(Max_t_d)
        Nst = len(Ntrace)
        GFmatrix = np.zeros((np.array(Ntrace, dtype=int).sum() + Max_t_d*Nst, 5))
        tb = 0
    else:
        GFmatrix = np.zeros((np.array(Ntrace, dtype=int).sum(), 5))

    #### Inversion:
    for i, trid in enumerate(trlist):
        trmeta = metadata[trid]
        trlat = trmeta['latitude']
        trlon = trmeta['longitude']

        try :
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
            synth = ltrim(synth, t_p - Max_t_d, delta)
            synth = synth[:, :Ntrace[i] + Max_t_d]
        else:
            synth = ltrim(synth, t_p - t_d, delta)
            synth = synth[:, :Ntrace[i]]

        if OnlyGetFullGF:
            ta = tb
            tb = ta + Max_t_d + (indexes[i+1] - indexes[i])
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
    inversion = lstsq(GFmatrix, ObservedDisp, rcond=None)
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
        misfit = 100.*np.sqrt(np.sum((syn-ObservedDisp)**2)/np.sum(syn**2))

    if GFs:
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



def remove_individual_traces(tol, M, GFmatrix, ObservedDisp,
                                  trlist, Ntrace):
    '''
    Remove any trace with an indivudual misfit higher than *tol*.

    :param float tol: The maximum acceptable value for a misfit.
    :param array M: Moment tensor.
    :param array GFmatrix: Greens function matrix.
    :param array Observed Disp: The observed displacements.
    :param list trlist: ???
    :param list Ntrace: List containing the lengths of the traces in *trlist*.

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
    Ntrace_fix = []
    M = np.delete(M,2)

    # construct the synthetic trace
    syn = GFmatrix.dot(M)
    j = 0
    for i, trid in enumerate(trlist):
        # get the observed and synthetic traces for a channel
        n1 = j
        n2 = Ntrace[i] + j
        obs = ObservedDisp[n1:n2]
        syn2 = syn[n1:n2]

        # calculate the misfit for the channel
        misfit = 100.*np.linalg.norm(syn2-obs)\
                 /np.linalg.norm(syn2)

        j += Ntrace[i]

        # add the trace to the remaining traces if misfit is within tolerance.
        if misfit <= tol :
            ObservedDisp_fix = np.append(ObservedDisp_fix, obs)
            GFmatrix_fix = np.append(GFmatrix_fix, GFmatrix[n1:n2,:], axis=0)
            trlist_fix.append(trid)
            Ntrace_fix.append(Ntrace[i])

    # build the result
    GFmatrix = GFmatrix_fix[1:,:]
    ObservedDisp = ObservedDisp_fix[1:]
    trlist = trlist_fix[:]
    Ntrace = Ntrace_fix[:]

    return GFmatrix, ObservedDisp, trlist, Ntrace


def rot_12_NE(st, META):
    '''
    Performs a  12 -> NE rotation.

    :param st: A stream containing the traces to rotate.
    :param dict META: Dictionary contaning the metadata for all streams.

    :return: The rotated streams.
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
            logger.warning("%s Channel 2 not found. Impossible to rotate", tr.id)
            continue

        timeA = max(tr1.stats.starttime,tr2.stats.starttime)
        timeB = min(tr1.stats.endtime,tr2.stats.endtime)
        tr1.trim(timeA,timeB)
        tr2.trim(timeA,timeB)
        azi = META[id1]['azimuth']
        tr1.data, tr2.data = Rot2D(tr1.data,tr2.data,-azi)

    return st






def Rot2D(x,y,angle):
    '''
    Given 2 vector x,y (same length) it performs a couterclockwise
    2d-rotation in a angle "angle" (degrees) for earch pair of elements
    '''
    ang = -angle*np.pi/180.
    xr =  x*np.cos(ang) - y*np.sin(ang)
    yr =  x*np.sin(ang) + y*np.cos(ang)
    return xr,yr



def get_timedelay_misfit_wrapper(args):
    return get_timedelay_misfit(*args)



def get_timedelay_misfit(t_d, GFmatrix, Ntrace, ObservedDisp, max_t_d):
    try:
        max_t_d = int(max_t_d)
        t_d = int(t_d)
        t_d2 = max_t_d - t_d
        GFmatrix_sm = np.zeros((np.array(Ntrace,dtype=np.int).sum(),5))
        cumNtraces = np.concatenate(([0], np.array(Ntrace,dtype=np.int).cumsum()))
        #aryNtrace = np.asarray(Ntrace, dtype=np.int)
        #GFmatrix_sm = np.zeros((aryNtrace.sum(),5))
        #cumNtraces = np.concatenate(([0], aryNtrace.cumsum()))
        tb_fm = 0
        for i_ntr, ta in enumerate(cumNtraces[:-1]):
            trlen = Ntrace[i_ntr]
            tb = ta + trlen
            ta_fm = tb_fm
            tb_fm = ta_fm +  max_t_d + tb-ta
            GFmatrix_sm[ta:tb] = GFmatrix[ta_fm + t_d2 :ta_fm + t_d2 + trlen]

        inversion = lstsq(GFmatrix_sm, ObservedDisp, rcond=None)
        return inversion[1][0]
    except Exception:
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))


def ltrim(data, starttime, delta):
    '''
    Left trimming similar to obspy but when *data* is a np array.

    :param data: array to trim, with time being the last axis.
    :type data: :py:class:`numpy.array`.
    :param numeric starttime: The amount of time to trim from the front of the trace in seconds.
    :param float delta: Sampling interval in seconds.

    :returns: A view of the relevant subset of *data*.
    '''
    i_of = int(round(starttime/delta))
    if i_of < 0:
        gap = data[..., 0, None] * np.ones(abs(i_of))
        return np.concatenate((gap, data), axis=-1)
    else:
        return data[..., i_of:]

def MT_result(M, misfit, depth, t_d):
    logger.info("Mrr: % e", M[0])
    logger.info("Mtt: % e", M[1])
    logger.info("Mpp: % e", M[2])
    logger.info("Mrt: % e", M[3])
    logger.info("Mrp: % e", M[4])
    logger.info("Mtp: % e", M[5])
    logger.info("misfit: %.0f%%", misfit)
    M2 = M*M
    m0 = np.sqrt(0.5 * (M2[0] + M2[1] + M2[2]) + M2[3] + M2[4] + M2[5])
    mag = 2./3.*(np.log10(m0)-9.10)
    logger.info("m0: % e", m0)
    logger.info("magnitude: %.5f", mag)
    return {
        'Mrr': M[0],
        'Mtt': M[1],
        'Mpp': M[2],
        'Mrt': M[3],
        'Mrp': M[4],
        'Mtp': M[5],
        'misfit': misfit,
        'm0': m0,
        'magnitude': round(mag, 1),
        'depth': depth,
        'time_delay': t_d,
    }
