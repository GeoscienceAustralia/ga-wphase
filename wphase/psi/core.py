'''
Core functionality for performing the Wphase inversion.
'''

import sys, os, glob, traceback
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
except:
    from obspy.core.util.geodetics import locations2degrees

try:
    from obspy.geodetics.base import gps2dist_azimuth
except:
    from obspy.core.util.geodetics import gps2DistAzimuth as gps2dist_azimuth

try:
    from obspy.signal.invsim import paz_2_amplitude_value_of_freq_resp
except:
    from obspy.signal.invsim import paz2AmpValueOfFreqResp as paz_2_amplitude_value_of_freq_resp

_greens_function_dir = None
_bpfunction = None
_N_TASKS_PER_PROC = None

def poolInitialiser(gfdir):
    """
    Initialise a worker in the pool. This is only required on windows where
    forking is not supported.
    """

    global _greens_function_dir
    _greens_function_dir = gfdir

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
    bpf=None,
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

    :param callable bpf: None or afunction to perform the Wphase filtering.
        This would typically be a buttworth filter. The signature must be as
        follows\:

        .. code-block:: python

            def bandpassfilter(sta, data, length, delta, corners, npass, low, high):
                pass

        where the arguments are\:

        - *sta*: The station name.
        - *data*: numpy array containing the data to be filtered.
        - *length*: The length of *data*.
        - *delta*: Sampling interval in seconds.
        - *corners*: Integer specifying the order of the filter.
        - *npass*: Integer specifying the number of passes: 1 for forward
            filtering only and 2 for forward and reverse (i.e. zero phase)
            filtering.
        - *low*: Low frequency cutoff of filter in hertz.
        - *high*: High frequency cutoff of filter in hertz.


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
            #. Inversion result for each grid point in inputs_latlon
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

    # set the function to use for filtering the Wphase
    global _bpfunction
    if bpf is None:
        from .bandpass import bandpassfilter
        def bpf(sta, data, length, delta, corners, npass, low, high):
            return bandpassfilter(data, delta, corners, low, high)
    _bpfunction = bpf

    # set the directory where the greens functions live
    global _greens_function_dir
    gfdir = os.path.normpath(gfdir)
    _greens_function_dir = gfdir #Directory of the Green Functions

    # check if we are using an hdf version of the greens functions
    if gfdir.split('.')[-1] == 'hdf5':
        hdf5_flag = True
        with h5py.File(gfdir , "r") as GFfile:
            hdirs_hdf5 = [dep for dep in GFfile]
    else:
        hdf5_flag = False
        hdirs_hdf5 = None

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

    print 'initial number of traces: {}'.format(len(st_sel))

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
        tr.data = np.array(tr.data, dtype=np.float)

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
            print "Unknown transfer function. Skipping", tr.id
            trlist.remove(tr.id)
            continue

        # Fitting the instrument response and getting coefficients
        (om0, h, G) = getCOEFFfit(trmeta['sensitivity'],freq,AmpfromPAZ)
        if not np.all(np.isfinite([om0, h, G])):
            print "Imposible to get Coeff. Skipping", tr.id
            trlist.remove(tr.id)
            continue

        AmpfromCOEFF= np.abs(omega*omega / \
                (omega*omega + 2j*h*om0*omega - om0*om0))

        # L2 norm:
        misfit = 100*np.linalg.norm(AmpfromPAZ-AmpfromCOEFF) \
                / np.linalg.norm(AmpfromPAZ)

        if misfit >  response_misfit_tol:
            print 'Bad fitting for response. Skipping:\n{}\t {: E}\n'.format(
                    tr.id,
                    misfit)
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
            print str(e)
            trlist.remove(tr.id)
            continue

        # trim to the Wphase time window
        tr.trim(t1,t2)
        if len(tr)== 0:
            print "Empty trace. Skipping", tr.id
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

    print "OL1:"
    print "Preliminary W-phase magnitude for the event:", np.round(pre_wp_mag, 1)

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
    if hdf5_flag:
        with h5py.File(_greens_function_dir , "r") as GFfile:
            dt = GFfile.attrs["dt"]
    else:
        dtDir = os.path.join(_greens_function_dir, "H003.5", "PP", "GF.0001.SY.LHZ.SAC")
        dt = read(dtDir)[0].stats.delta

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

        tr.data = np.array(tr.data,dtype=np.float)
        if trmeta['transfer_function'] == "B":
            AmpfromPAZ  = Vpaz2freq(trmeta,freq/2./np.pi)  # hz to rad*hz
        elif trmeta['transfer_function'] == "A":
            AmpfromPAZ  = Vpaz2freq(trmeta,freq)
        else:
            print "Unknown transfer function. Skipping", tr.id
            trlist.remove(tr.id)
            continue

        (om0, h, G) = getCOEFFfit(trmeta['sensitivity'],freq,AmpfromPAZ)
        if not np.all(np.isfinite([om0, h, G])):
            print "Imposible to get Coeff. Skipping", tr.id
            trlist.remove(tr.id)
            continue

        DATA_INFO[tr.id] = [om0, h, G, dt, t1, t2]

        try:
            tr.data, coeff = RTdeconv(
                tr,
                om0,
                h,
                G,
                dt,corners=4,
                baselinelen=60./dt,
                taperlen= 10.,
                fmin = 1./Tb,
                fmax = 1./Ta,
                get_coef = True)

        except RTdeconvError as e:
            print str(e)
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

    pool=Pool(processes, poolInitialiser, (_greens_function_dir,), _N_TASKS_PER_PROC)

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
        hdf5_flag,
        OnlyGetFullGF=True,
        Max_t_d=max_t_d,
        hdirs_hdf5=hdirs_hdf5)

    # inputs for the TIME DELAY search
    inputs = [(t_d_test, GFmatrix, Ntrace, ObservedDisp, max_t_d) for t_d_test in time_delays]
    misfits = pool.map(get_timedelay_misfit_wrapper,inputs)
    pool.close()
    pool.join()

    # Set t_d (time delay) and and t_h (half duration) to optimal values:
    mis_min =  np.array(misfits).argmin()
    output_dic['misfits'] = {'min':mis_min, 'array':misfits}
    t_d = t_h = time_delays[mis_min]
    MRF = MomentRateFunction(t_h, dt)
    print "Source time function, time delay:", len(MRF), t_d

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
        hdf5_flag,
        GFs=True,
        hdirs_hdf5=hdirs_hdf5)

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
            print msg
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
            hdf5_flag,
            GFs=True,
            hdirs_hdf5=hdirs_hdf5)

    syn = (M[0]*GFmatrix[:,0] + M[1]*GFmatrix[:,1] +
           M[3]*GFmatrix[:,2] + M[4]*GFmatrix[:,3] +
           M[5]*GFmatrix[:,4])

    print "OL2:"
    print "Mrr: ", '{: e}'.format(M[0])
    print "Mtt: ", '{: e}'.format(M[1])
    print "Mpp: ", '{: e}'.format(M[2])
    print "Mrt: ", '{: e}'.format(M[3])
    print "Mrp: ", '{: e}'.format(M[4])
    print "Mtp: ", '{: e}'.format(M[5])

    print "misfit: ", misfit

    M2 = M*M
    m0 = np.sqrt(0.5 * (M2[0] + M2[1] + M2[2]) + M2[3] + M2[4] + M2[5])
    mag = 2./3.*(np.log10(m0)-9.10)
    print "m0: ", m0
    print "magnitude: ", mag

    output_dic['OL2'] = {}
    output_dic['OL2']['Mrr']= M[0]
    output_dic['OL2']['Mtt']= M[1]
    output_dic['OL2']['Mpp']= M[2]
    output_dic['OL2']['Mrt']= M[3]
    output_dic['OL2']['Mrp']= M[4]
    output_dic['OL2']['Mtp']= M[5]
    output_dic['OL2']['misfit']= misfit
    output_dic['OL2']['m0']= m0
    output_dic['OL2']['magnitude']= round(mag,1)
    output_dic['OL2']['depth'] = hypdep
    output_dic['OL2']['time_delay'] = t_d

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

    #~ deps_grid = get_depths_for_grid(hypdep)
    #~ lat_grid, lon_grid = get_latlon_for_grid(hyplat,hyplon)
    #~
    #~ deps_grid = get_depths_for_grid(hypdep)
    #~ lat_grid, lon_grid = get_latlon_for_grid(hyplat,hyplon,dist_lat = 3.0,
                            #~ dist_lon = 3.0,delta = 0.4)
    ##Testing shorcut:
    ##latlon search
    lat_grid, lon_grid = get_latlon_for_grid(hyplat,hyplon,dist_lat = 3.0,
                                             dist_lon = 3.0,delta = 0.8)
    #lat_grid = lat_grid[::2]
    #lon_grid = lon_grid[::2]
    pool=Pool(processes, poolInitialiser, (_greens_function_dir,), _N_TASKS_PER_PROC)
    inputs_latlon = [(t_h,t_d,(hyplatg,hyplong,hypdep), orig, (Ta,Tb),
              MRF, ObservedDisp, Ntrace, metadata, trlist, hdf5_flag,
              {'residuals':False,'hdirs_hdf5':hdirs_hdf5})
              for hyplatg in lat_grid
              for hyplong in lon_grid]
    latlon_search = pool.map(core_inversion_wrapper,inputs_latlon)
    pool.close()
    pool.join()

    #latlon_search = map(core_inversion_wrapper,inputs_latlon)
    misfits_latlon = np.array([latlon_search[i][1] for i in range(len(inputs_latlon))])
    cenlat, cenlon = inputs_latlon[misfits_latlon.argsort()[0]][2][:2]
    moments = latlon_search  #Compatibility

    if hdf5_flag:
        deps_grid = get_depths_for_grid_hdf5(hypdep)
    else:
        deps_grid = get_depths_for_grid(hypdep)#[::2]

    pool=Pool(processes, poolInitialiser, (_greens_function_dir,), _N_TASKS_PER_PROC)
    inputs_dep = [(t_h,t_d,(cenlat,cenlon,hypdepg), orig, (Ta,Tb),
              MRF, ObservedDisp, Ntrace, metadata, trlist, hdf5_flag,
              {'residuals':False,'hdirs_hdf5':hdirs_hdf5})
              for hypdepg in deps_grid]

    depth_search = pool.map(core_inversion_wrapper,inputs_dep)
    pool.close()
    pool.join()

    #depth_search = map(core_inversion_wrapper,inputs_dep)
    misfits_depth = np.array([depth_search[i][1] for i in range(len(inputs_dep))])
    cendep = deps_grid[misfits_depth.argsort()[0]]


    #########old grid search
    #~ deps_grid = get_depths_for_grid(hypdep)[::2]
    #~ lat_grid, lon_grid = get_latlon_for_grid(hyplat,hyplon,dist_lat = 3.0,
                            #~ dist_lon = 3.0,delta = 0.4)
    #~ lat_grid = lat_grid[::2]
    #~ lon_grid = lon_grid[::2]
    #~
    #~ pool=Pool(processes, poolInitialiser, (_greens_function_dir,), _N_TASKS_PER_PROC)
    #~ inputs = [(t_h,t_d,(hyplatg,hyplong,hypdepg), orig, (Ta,Tb),
              #~ MRF, ObservedDisp, Ntrace, metadata, trlist, {'residuals':False})
              #~ for hypdepg in deps_grid
              #~ for hyplatg in lat_grid
              #~ for hyplong in lon_grid]
            #~
    #~ moments = pool.map(core_inversion_wrapper,inputs)
    #~ misfits = [moments[i][1] for i in range(len(inputs))]
    #~ min_misfit_args = np.array(misfits).argsort()[0:Nfinersearch]
    #~ min_misfit_arg = np.array(misfits).argmin()
    #~ ####Second (finer) gridsearch:
    #~ global_min_misfit  = 10. # Initial Value, should be > 1.
    #~ ## If Nfinersearch = 0 I will not do a finner grid search
    #~ if len(min_misfit_args) == 0:
        #~ cenlat = inputs[min_misfit_arg][2][0]
        #~ cenlon = inputs[min_misfit_arg][2][1]
        #~ cendep = inputs[min_misfit_arg][2][2]
    #~ for min_misfit_arg in min_misfit_args:
        #~ #### Redefining the hyp loc:
        #~ hyplat = inputs[min_misfit_arg][2][0]
        #~ hyplon = inputs[min_misfit_arg][2][1]
        #~ hypdep = inputs[min_misfit_arg][2][2]
        #~ print "Looking for an optimal centroid location close to:"
        #~ print hyplat,hyplon,hypdep
        #~
        #~ lat_grid, lon_grid = get_latlon_for_grid(hyplat,hyplon, dist_lat = 0.3,
                                #~ dist_lon = 0.3,delta = 0.1 )
        #~ deps_grid = get_depths_for_grid(hypdep, length = 30)
        #~ pool=Pool(processes, poolInitialiser, (_greens_function_dir,), _N_TASKS_PER_PROC)
        #~ inputs_f = [(t_h,t_d,(hyplatg,hyplong,hypdepg), orig, (Ta,Tb),
                #~ MRF, ObservedDisp, Ntrace, metadata, trlist, {'residuals':False})
                #~ for hypdepg in deps_grid
                #~ for hyplatg in lat_grid
                #~ for hyplong in lon_grid]
                #~
        #~ moments = pool.map(core_inversion_wrapper,inputs_f)
        #~ misfits = [moments[i][1] for i in range(len(inputs_f))]
        #~ misfits = np.array(misfits)
        #~ if misfits.min() < global_min_misfit:
            #~ min_misfit_arg = misfits.argmin()
            #~ cenlat = inputs_f[min_misfit_arg][2][0]
            #~ cenlon = inputs_f[min_misfit_arg][2][1]
            #~ cendep = inputs_f[min_misfit_arg][2][2]
            #~ print "Found:", cenlat, cenlon, cendep
            #~ global_min_misfit = misfits.min()
        #~ else:
            #~ cenlat = inputs[min_misfit_arg][2][0]
            #~ cenlon = inputs[min_misfit_arg][2][1]
            #~ cendep = inputs[min_misfit_arg][2][2]
    ###############End old grid search

    ###Final inversion!!

    M, misfit, GFmatrix = core_inversion(t_h,t_d,(cenlat,cenlon,cendep),
                                orig, (Ta,Tb), MRF, ObservedDisp,
                                 Ntrace, metadata, trlist, hdf5_flag,
                                 GFs = True,hdirs_hdf5=hdirs_hdf5)

    syn = (M[0]*GFmatrix[:,0] + M[1]*GFmatrix[:,1]
          + M[3]*GFmatrix[:,2] + M[4]*GFmatrix[:,3]
          + M[5]*GFmatrix[:,4])


    print "OL3:"
    print "Mrr: ", '{: e}'.format(M[0])
    print "Mtt: ", '{: e}'.format(M[1])
    print "Mpp: ", '{: e}'.format(M[2])
    print "Mrt: ", '{: e}'.format(M[3])
    print "Mrp: ", '{: e}'.format(M[4])
    print "Mtp: ", '{: e}'.format(M[5])

    print "misfit: ", misfit

    M2 = M*M
    m0 = np.sqrt(0.5*(M2[0]+M2[1]+M2[2])+M2[3]+M2[4]+M2[5])
    mag = 2./3.*(np.log10(m0)-9.10)
    print "m0: ", m0
    print "magnitude: ", mag
    cenloc = (cenlat,cenlon,cendep)

    # QZ
    output_dic['OL3'] = {}
    output_dic['OL3']['Mrr']= M[0]
    output_dic['OL3']['Mtt']= M[1]
    output_dic['OL3']['Mpp']= M[2]
    output_dic['OL3']['Mrt']= M[3]
    output_dic['OL3']['Mrp']= M[4]
    output_dic['OL3']['Mtp']= M[5]
    output_dic['OL3']['misfit']= misfit
    output_dic['OL3']['m0']= m0
    output_dic['OL3']['magnitude']= round(mag,1)
    output_dic['OL3']['depth'] = cendep
    output_dic['OL3']['time_delay'] = t_d

    '''
    print "****************   output_dic   *****************"
    print output_dic
    print "****************   output_dic   *****************"
    '''

    return M, ObservedDisp, syn, trlist, Ntrace, cenloc, t_d, inputs_latlon, moments, DATA_INFO


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

    accel = _bpfunction(sta,accel,len(accel),dt ,corners , 1 , fmin,fmax)

    vel = np.zeros(len(data))
    vel[1:] = 0.5*dt*np.cumsum(accel[:-1]+accel[1:])

    dis = np.zeros(len(data))
    dis[1:] = 0.5*dt*np.cumsum(vel[:-1]+vel[1:])

    if get_coef:
        return dis, (c0,c1,c2)

    return dis



def GFSelectZ(dist, hypdepth, GFdir, *args, **kwargs):
    """
    Loads the Greens function for the given *dist* and *hypdepth* from a SAC file.

    :param float dist: The epicentral distance.
    :param float hypdepth: The depth.

    :return: Tuple containin the four greens functions corresponding to the vertical component.
    :rtype: Tuple of four arrays.
    """

    ###I'm considering the decimal .5 in the depths dirs.

    depth = []
    depthdir = []
    for file in os.listdir(GFdir):
        if file[-2:] == ".5":
            depthdir.append(file)
            depth.append(float(file[1:]))
    BestDirIndex = np.argsort(abs(hypdepth-np.array(depth)))[0]
    # hdir is the absolute path to the closest deepth.
    hdir = os.path.join(GFdir, depthdir[BestDirIndex])
    dist_dir =  int(dist*10.+0.5)

    if dist_dir > 900:
        dist_dir = 900 # We do not have GFs for dist > 90 deg
    dist_str = str(dist_dir)#Some GFs have only odd dists.


    #dist_str = str(int(dist*10.))#Some GFs have only odd dists.
    dist_form = dist_str.zfill(4)

    ## Loading files
    trPPPath = os.path.join(hdir, "PP", "GF." + dist_form + ".SY.LHZ.SAC")
    trPP = read(trPPPath)[0]
    #trPP.data -= trPP.data[0]

    trRRPath = os.path.join(hdir, "RR", "GF." + dist_form + ".SY.LHZ.SAC")
    trRR = read(trRRPath)[0]
    #trRR.data -= trRR.data[0]


    trRTPath = os.path.join(hdir, "RT", "GF." + dist_form + ".SY.LHZ.SAC")
    trRT = read(trRTPath)[0]
    #trRT.data -= trRT.data[0]

    trTTPath = os.path.join(hdir, "TT", "GF." + dist_form + ".SY.LHZ.SAC")
    trTT = read(trTTPath)[0]
    #trTT.data -= trTT.data[0]

    return (trPP, trRR, trRT,  trTT)



def MTrotationZ(azi, (trPPsy,  trRRsy, trRTsy,  trTTsy)):
    """
    Rotates the Greens function for the given azimuth (*azi*).

    :param float azi: The station azimuth with respect to the hypocenter.
    :param (trPPsy,  trRRsy, trRTsy,  trTTsy): Output of :py:func:`GFSelectZ`.

    :return: Tuple containin the four greens functions corresponding to the vertical component.
    :rtype: Tuple of six arrays.
    """

    sinp = np.sin(azi)
    cosp = np.cos(azi)


	#Creating traces for the rotated syn:
    trRRsy_rot =   trRRsy.copy()
    trPPsy_rot =   trPPsy.copy()
    trTTsy_rot =   trTTsy.copy()
    trTPsy_rot =   trTTsy.copy()
    trRTsy_rot =   trRTsy.copy()
    trRPsy_rot =   trPPsy.copy()

    #Obtaining the rotated synthetics:
    trRRsy_rot.data =   trRRsy.data
    trPPsy_rot.data =   sinp*sinp*trTTsy.data+cosp*cosp*trPPsy.data
    trTTsy_rot.data =   cosp*cosp*trTTsy.data+sinp*sinp*trPPsy.data
    trTPsy_rot.data =   2.*sinp*cosp*(trTTsy.data-trPPsy.data)
    trRTsy_rot.data =   cosp*trRTsy.data
    trRPsy_rot.data =   sinp*trRTsy.data

    return [trRRsy_rot, trPPsy_rot,  trTTsy_rot, trTPsy_rot,trRTsy_rot, trRPsy_rot]


def core_inversion(t_h,t_d,cmtloc,orig, periods, MRF,
                      ObservedDisp, Ntrace, metadata, trlist,
                      hdf5_flag, GFs = False, residuals = False,
                      OnlyGetFullGF=False, Max_t_d = 200, hdirs_hdf5= None):
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
    :param bool hdf5_flag: *True* if the greens functions are stored in HDF5, *False* if they are in SAC.
    :param bool GFs: *True* if the greens functions should be returned.
    :param bool residuals: *True*, return the 'raw' misfit from the least squares inversion,
        *False*, return the 'relative' misfit as a percentage: i.e. the 'raw' misfit divided by
        the norm of the synthetics.
    :param bool OnlyGetFullGF: *True*, return the greens function matrix for the maximum time delay (*Max_t_d*)
        without performing the inversion, *False*, perform the inversion for the time delay given
        by *t_d*.
    :param numeric Max_t_d: Maximum time delay to consider if *OnlyGetFullGF = True*.
    :param hdirs_hdf5: List of the names of the depth directories in the greens function HDF5 library.
        Ignored if *hdf5_flag = False*.

    :return: What is returned depends on the values of the parameters as described below.

        - If *OnlyGetFullGF = True*, then just the greens function matrix.
        - Otherwise a tuple containing:

            0. moment tensor components in Nm ['RR', 'PP', 'TT', 'TP', 'RT', 'RP']
            #. misfit (percent L2 norm misfit error of the solution), and

            If *GFs = True*

            2. The greens function matrix also.
    '''

    # Dealing with hdf5 files -  define the functions used to extract greens functions.
    if hdf5_flag:
        GFSelectZ_ = GFSelectZ_hdf5
        MTrotationZ_ = MTrotationZ_hdf5
        GFSelectN_ = GFSelectN_hdf5
        MTrotationN_ = MTrotationN_hdf5
        GFSelectE_ = GFSelectE_hdf5
        MTrotationE_ = MTrotationE_hdf5
        with h5py.File(_greens_function_dir , "r") as GFfile:
            delta = GFfile.attrs["dt"]
    else:
        GFSelectZ_ = GFSelectZ
        MTrotationZ_ = MTrotationZ
        GFSelectN_ = GFSelectN
        MTrotationN_ = MTrotationN
        GFSelectE_ = GFSelectE
        MTrotationE_ = MTrotationE
        dtDir = os.path.join(_greens_function_dir, "H003.5", "PP", "GF.0001.SY.LHZ.SAC")
        delta = read(dtDir)[0].stats.delta

    hyplat = cmtloc[0]
    hyplon = cmtloc[1]
    hypdep = cmtloc[2]
    Ta     = periods[0]
    Tb     = periods[1]

    # Index of the first value that will be valid after convolution with MRF.
    FirstValid = int(len(MRF)/2.)

    # the indicies of the beginning and end of each trace in observed displacements (ObservedDisp)
    indexes =  np.array(np.concatenate((np.array([0.]), np.cumsum(Ntrace))), dtype='int')


    if OnlyGetFullGF:
        Max_t_d = int(Max_t_d)
        Nst = len(Ntrace)
        GFmatrix = np.zeros((np.array(Ntrace, dtype=np.int).sum() + Max_t_d*Nst, 5))
        tb = 0
    else:
        GFmatrix = np.zeros((np.array(Ntrace, dtype=np.int).sum(), 5))

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
        if trid[-1] == 'Z':
            syntr = GFSelectZ_(dist, hypdep, _greens_function_dir, hdirs_hdf5)
            syntrrot = MTrotationZ_(-azi*np.pi/180, syntr)
        else:
            syntrT = GFSelectN_(dist, hypdep, _greens_function_dir, hdirs_hdf5)
            syntrP = GFSelectE_(dist, hypdep, _greens_function_dir, hdirs_hdf5)
            syntrrotT = MTrotationN_(-azi*np.pi/180, syntrT)
            syntrrotP = MTrotationE_(-azi*np.pi/180, syntrP)

            for l in range(6):
                if hdf5_flag:
                    syntrrotT[l], syntrrotP[l] = Rot2D(-syntrrotT[l], syntrrotP[l], -bazi)
                else:
                    syntrrotT[l].data, syntrrotP[l].data = Rot2D(-syntrrotT[l].data,
                                                            syntrrotP[l].data, -bazi)
            if trid[-1] == '1' or trid[-1] == 'N':
                syntrrot = syntrrotT


            elif trid[-1] == '2' or trid[-1] == 'E':
                syntrrot = syntrrotP

        # Extract Wphase from synthetics
        for i_trR, trR in enumerate(syntrrot): ##Check this loop, seems to be slow.
            if hdf5_flag:
                data_ = trR
            else:
                data_ = trR.data

            # greens functions are not stored in standard units, convert ot Nm (Newton meters)
            data_ *= 10.**-21

            # use only what we think is noise to calculate the mean
            mean = np.mean(data_[:60])
            data_ -= mean

            data_ = ndimage.convolve(data_, MRF, mode='nearest')\
                            [FirstValid-2:-FirstValid]

            # TODO: look at this... do we want to remove this catch?
            try:
                fillvec1 = data_[ 0] * np.ones(FirstValid)
                fillvec2 = data_[-1] * np.ones(FirstValid)
            except:
                print FirstValid

            data_ = np.concatenate((fillvec1, data_, fillvec2))
            mean = np.mean(data_[:60])
            data_ -= mean
            data_ = _bpfunction(sta,data_,len(trR),delta ,4 ,1 , 1./Tb, 1./Ta)

            if OnlyGetFullGF:
                data_ = ltrim(data_, t_p - Max_t_d, delta)
                data_ = data_[:Ntrace[i] + Max_t_d ]
            else:
                data_ = ltrim(data_, t_p - t_d, delta)
                data_ = data_[:Ntrace[i]]
            if hdf5_flag:
                syntrrot[i_trR]  = data_
            else:
                trR.data = data_

        if OnlyGetFullGF:
            ta = tb
            tb = ta + Max_t_d + (indexes[i+1] - indexes[i])
        else:
            ta = indexes[i]
            tb = indexes[i+1]

        # the first of the following lines are due to constraint that volume does not change
        GFmatrix[ta:tb, 0] = syntrrot[0][:] - syntrrot[1][:]
        GFmatrix[ta:tb, 1] = syntrrot[2][:] - syntrrot[1][:]
        GFmatrix[ta:tb, 2] = syntrrot[4][:]
        GFmatrix[ta:tb, 3] = syntrrot[5][:]
        GFmatrix[ta:tb, 4] = syntrrot[3][:]

    if OnlyGetFullGF:
        return GFmatrix

    # perform the inversion
    inversion = lstsq(GFmatrix, ObservedDisp, rcond=None)
    M = inversion[0]

    # construct the synthetics
    syn = (M[0]*GFmatrix[:,0] + M[1]*GFmatrix[:,1]
        + M[2]*GFmatrix[:,2] + M[3]*GFmatrix[:,3] + M[4]*GFmatrix[:,4])

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



def get_depths_for_grid(hypdep, length = 100):
    '''
    This function will return a list with the adequate values
    of depth to look at in the grid search
    '''

    ####Parameters:
    min_dep = 11. #No shallower depths allowed.
     # Total length of the section in which we'll search
    #####
    depths = sorted([float(dir.split(os.path.sep)[-1][1:])
              for dir in glob.glob(os.path.join(_greens_function_dir, '*.5'))])
    depths = np.array(depths)
    depths_grid = depths[np.where(np.abs(depths-hypdep)<length*0.5)]
    depths_grid = depths_grid[np.where(depths_grid > min_dep)]
    return depths_grid



def get_depths_for_grid_hdf5(hypdep, length = 100):
    '''
    This function will return a list with the adequate values
    of depth to look at in the grid search
    '''

    ####Parameters:
    min_dep = 11. #No shallower depths allowed.
     # Total length of the section in which we'll search
    #####
    with h5py.File(_greens_function_dir , "r") as GFfile:
        depths = sorted(GFfile.keys())
        depths = np.array([dep[1:] for dep in depths], dtype=np.float)
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
    except:
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
    syn = (M[0]*GFmatrix[:,0] + M[1]*GFmatrix[:,1]+ M[2]*GFmatrix[:,2]
             + M[3]*GFmatrix[:,3] + M[4]*GFmatrix[:,4])
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
            print tr.id, "Channel 2 not found. Impossible to rotate"
            continue

        timeA = max(tr1.stats.starttime,tr2.stats.starttime)
        timeB = min(tr1.stats.endtime,tr2.stats.endtime)
        tr1.trim(timeA,timeB)
        tr2.trim(timeA,timeB)
        azi = META[id1]['azimuth']
        tr1.data, tr2.data = Rot2D(tr1.data,tr2.data,-azi)

    return st



def GFSelectN(dist, hypdepth, GFdir, *args, **kwargs):
    ###I'm considering the decimal .5 in the depths dirs.
    depth = []
    depthdir = []
    for file in os.listdir(GFdir):
        if file[-2:] == ".5":
            depthdir.append(file)
            depth.append(float(file[1:]))
    BestDirIndex = np.argsort(abs(hypdepth-np.array(depth)))[0]
    # hdir is the absolute path to the closest deepth.
    hdir = os.path.join(GFdir, depthdir[BestDirIndex])
    dist_dir =  int(dist*10.+0.5)

    if dist_dir > 900:
        dist_dir = 900 # We do not have GFs for dist > 90 deg
    dist_str = str(dist_dir)#Some GFs have only odd dists.
    dist_form = dist_str.zfill(4)

    ## Loading files
    trPP = read(os.path.join(hdir, "PP", "GF." + dist_form + ".SY.LHL.SAC"))[0]
    trRR = read(os.path.join(hdir, "RR", "GF." + dist_form + ".SY.LHL.SAC"))[0]
    trRT = read(os.path.join(hdir, "RT", "GF." + dist_form + ".SY.LHL.SAC"))[0]
    trTT = read(os.path.join(hdir, "TT", "GF." + dist_form + ".SY.LHL.SAC"))[0]

    return (trPP, trRR, trRT,  trTT)



def MTrotationN(azi, (trPPsy,  trRRsy, trRTsy,  trTTsy)):
    '''
    Rotate the Moment Tensor according to the azimuthal angle in radians.
    '''

    sinp = np.sin(azi)
    cosp = np.cos(azi)

	#Creating traces for the rotated syn:
    trRRsy_rot =   trRRsy.copy()
    trPPsy_rot =   trPPsy.copy()
    trTTsy_rot =   trTTsy.copy()
    trTPsy_rot =   trTTsy.copy()
    trRTsy_rot =   trRTsy.copy()
    trRPsy_rot =   trPPsy.copy()

    #Obtaining the rotated synthetics:
    trRRsy_rot.data =   trRRsy.data
    trPPsy_rot.data =   sinp*sinp*trTTsy.data+cosp*cosp*trPPsy.data
    trTTsy_rot.data =   cosp*cosp*trTTsy.data+sinp*sinp*trPPsy.data
    trTPsy_rot.data =   2.*sinp*cosp*(trTTsy.data-trPPsy.data)
    trRTsy_rot.data =   cosp*trRTsy.data
    trRPsy_rot.data =   sinp*trRTsy.data

    return [trRRsy_rot, trPPsy_rot,  trTTsy_rot, trTPsy_rot,trRTsy_rot, trRPsy_rot]



def GFSelectE(dist, hypdepth, GFdir, *args, **kwargs):
    ###I'm considering the decimal .5 in the depths dirs.
    depth = []
    depthdir = []
    for file in os.listdir(GFdir):
        if file[-2:] == ".5":
            depthdir.append(file)
            depth.append(float(file[1:]))
    BestDirIndex = np.argsort(abs(hypdepth-np.array(depth)))[0]
    # hdir is the absolute path to the closest depth.
    hdir = os.path.join(GFdir, depthdir[BestDirIndex])
    dist_dir =  int(dist*10.+0.5)

    if dist_dir > 900:
        dist_dir = 900 # We do not have GFs for dist > 90 deg
    dist_str = str(dist_dir)#Some GFs have only odd dists.
    dist_form = dist_str.zfill(4)

    ## Loading files
    trRP = read(os.path.join(hdir, "RP", "GF." + dist_form + ".SY.LHT.SAC"))[0]
    trTP = read(os.path.join(hdir, "TP", "GF." + dist_form + ".SY.LHT.SAC"))[0]
    #trRP.data = -trRP.data
    #trTP.data = -trTP.data

    return (trRP, trTP)



def MTrotationE(azi, (trRPsy, trTPsy)):
    '''
    Rotates the Moment Tensor according to the azimuthal angle in radians.
    '''

    sinp = np.sin(azi)
    cosp = np.cos(azi)
    sin2p = np.sin(2.*azi)
    cos2p = np.cos(2.*azi)
	#Creating traces for the rotated syn:
    trRRsy_rot =   trRPsy.copy()
    trPPsy_rot =   trRPsy.copy()
    trTTsy_rot =   trRPsy.copy()
    trTPsy_rot =   trRPsy.copy()
    trRTsy_rot =   trRPsy.copy()
    trRPsy_rot =   trRPsy.copy()
    #Obtaining the rotated synthetics:
    trRRsy_rot.data[:] =   0.
    trPPsy_rot.data =   0.5*sin2p*trTPsy.data
    trTTsy_rot.data =  -0.5*sin2p*trTPsy.data
    trTPsy_rot.data =   cos2p*trTPsy.data
    trRTsy_rot.data =   -sinp*trRPsy.data
    trRPsy_rot.data =  cosp*trRPsy.data

    return [trRRsy_rot, trPPsy_rot,  trTTsy_rot, trTPsy_rot,trRTsy_rot, trRPsy_rot]



def MTrot(azi,mt):
    azi = np.pi/180.* azi
    sinp = np.sin(azi)
    cosp = np.cos(azi)
    sin2p = np.sin(2.*azi)
    cos2p = np.cos(2.*azi)

    mtrotpp = mt[1]*cosp**2-2*mt[3]*sinp*cosp+mt[2]*sinp**2
    mtrottt = mt[1]*sinp**2+2*mt[3]*sinp*cosp+mt[2]*cosp**2
    mtrotrr = float(mt[0])
    mtrottp = -0.5*(mt[2]-mt[1])*sin2p+mt[3]*cos2p
    mtrotrp = mt[5]*cosp - mt[4]*sinp
    mtrotrt = mt[4]*cosp + mt[5]*sinp

    return [mtrotrr, mtrotpp, mtrottt, mtrottp, mtrottp, mtrotrt, mtrotrp]



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
        tb_fm = 0
        for i_ntr, ta in enumerate(cumNtraces[:-1]):
            trlen = Ntrace[i_ntr]
            tb = ta + trlen
            ta_fm = tb_fm
            tb_fm = ta_fm +  max_t_d + tb-ta
            GFmatrix_sm[ta:tb] = GFmatrix[ta_fm + t_d2 :ta_fm + t_d2 + trlen]

        inversion = lstsq(GFmatrix_sm, ObservedDisp, rcond=None)
        return inversion[1][0]
    except:
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))



def GFSelectZ_hdf5(dist, hdir, gfdir, hdirs_hdf5):
    """
    Loads the Greens function for the given *dist* and *hypdepth* from a HDF file.

    :param float dist: The epicentral distance.
    :param float hypdepth: The depth.
    :param str gfdir: T HDF5 file.
    :param hdirs_hdf5: List of the names of the depth 'directories' inside the HDF5 file.
    :type hdirs_hdf5: List of str

    :return: Tuple containin the four greens functions corresponding to the vertical component.
    :rtype: Tuple of four arrays.
    """

    hdir = search_best_hdir_hdf5(hdir,hdirs_hdf5)#This can be done in core inv!
    #We do not have GF> 90. This is only to prevent the code from crashing!
    if dist > 90:
        dist = 90
    ###############################################################
    #Updated to work with the 0.1d spaced GF database
    #dist_str = str(int(dist*10.)/2*2+1)
    dist_str =  str(int(dist*10.+0.5))
    dist_form = dist_str.zfill(4)
    with h5py.File(gfdir , "r") as GFfile:
        trPP = GFfile[hdir + "PP/GF." + dist_form + ".SY.LHZ.SAC"][...]
        trRR = GFfile[hdir + "RR/GF." + dist_form + ".SY.LHZ.SAC"][...]
        trRT = GFfile[hdir + "RT/GF." + dist_form + ".SY.LHZ.SAC"][...]
        trTT = GFfile[hdir + "TT/GF." + dist_form + ".SY.LHZ.SAC"][...]

    return (trPP, trRR, trRT,  trTT)



def MTrotationZ_hdf5(azi, (trPPsy,  trRRsy, trRTsy,  trTTsy)):
    """
    Rotates the Greens function for the given azimuth (*azi*).

    :param float azi: The station azimuth with respect to the hypocenter.
    :param (trPPsy,  trRRsy, trRTsy,  trTTsy): Output of :py:func:`GFSelectZ`.

    :return: Tuple containin the four greens functions corresponding to the vertical component.
    :rtype: Tuple of six arrays.
    """

    #Obtaining the rotated synthetics:
    sinp = np.sin(azi)
    cosp = np.cos(azi)

    trRRsy_rot =   trRRsy
    trPPsy_rot =   sinp*sinp*trTTsy + cosp*cosp*trPPsy
    trTTsy_rot =   cosp*cosp*trTTsy+sinp*sinp*trPPsy
    trTPsy_rot =   2.*sinp*cosp*(trTTsy-trPPsy)
    trRTsy_rot =   cosp*trRTsy
    trRPsy_rot =   sinp*trRTsy
    return [trRRsy_rot, trPPsy_rot,  trTTsy_rot,
            trTPsy_rot,trRTsy_rot, trRPsy_rot]



def GFSelectN_hdf5(dist, hdir, gfdir, hdirs_hdf5):
    """
    Loads the Greens function for the given *dist* from a HDF file.

    :param float dist: The epicentral distance.
    :param str gfdir: T HDF5 file.
    :param hdirs_hdf5: List of the names of the depth 'directories' inside the HDF5 file.
    :type hdirs_hdf5: List of str

    :return: Tuple containin the four greens functions corresponding to the vertical component.
    :rtype: Tuple of four arrays.
    """

    ###I'm considering the decimal .5 in the depths dirs.
    hdir = search_best_hdir_hdf5(hdir,hdirs_hdf5)#This can be done in core inv!
    #We do not have GF> 90. This is only to prevent the code from crashing!
    if dist > 90:
        dist = 90
    ###############################################################
    dist_str =  str(int(dist*10.+0.5))
    dist_form = dist_str.zfill(4)

    ## Loading files
    with h5py.File(gfdir , "r") as GFfile:
        trPP = GFfile[hdir + "PP/GF." + dist_form + ".SY.LHL.SAC" ][...]
        trRR = GFfile[hdir + "RR/GF." + dist_form + ".SY.LHL.SAC" ][...]
        trRT = GFfile[hdir + "RT/GF." + dist_form + ".SY.LHL.SAC" ][...]
        trTT = GFfile[hdir + "TT/GF." + dist_form + ".SY.LHL.SAC" ][...]
    return (trPP, trRR, trRT,  trTT)



def MTrotationN_hdf5(azi, (trPPsy,  trRRsy, trRTsy,  trTTsy)):
    """
    Rotates the Greens function for the given azimuth (*azi*).

    :param float azi: The station azimuth with respect to the hypocenter.
    :param (trPPsy,  trRRsy, trRTsy,  trTTsy): Output of :py:func:`GFSelectZ`.

    :return: Tuple containin the four greens functions corresponding to the vertical component.
    :rtype: Tuple of six arrays.
    """

    sinp = np.sin(azi)
    cosp = np.cos(azi)
    #Obtaining the rotated synthetics:
    trRRsy_rot =   trRRsy
    trPPsy_rot =   sinp*sinp*trTTsy+cosp*cosp*trPPsy
    trTTsy_rot =   cosp*cosp*trTTsy+sinp*sinp*trPPsy
    trTPsy_rot =   2.*sinp*cosp*(trTTsy-trPPsy)
    trRTsy_rot =   cosp*trRTsy
    trRPsy_rot =   sinp*trRTsy

    return [trRRsy_rot, trPPsy_rot,  trTTsy_rot, trTPsy_rot,trRTsy_rot, trRPsy_rot]



def GFSelectE_hdf5(dist, hdir,gfdir,hdirs_hdf5):
    """
    Loads the Greens function for the given *dist* from a HDF file.

    :param float dist: The epicentral distance.
    :param str gfdir: T HDF5 file.
    :param hdirs_hdf5: List of the names of the depth 'directories' inside the HDF5 file.
    :type hdirs_hdf5: List of str

    :return: Tuple containin the four greens functions corresponding to the vertical component.
    :rtype: Tuple of four arrays.
    """

    ###I'm considering the decimal .5 in the depths dirs.
    hdir = search_best_hdir_hdf5(hdir,hdirs_hdf5)#This can be done in core inv!
    #We do not have GF> 90. This is only to prevent the code from crashing!
    if dist > 90:
        dist = 90
    ###############################################################
    dist_str =  str(int(dist*10.+0.5))
    dist_form = dist_str.zfill(4)
    ## Loading files
    with h5py.File(gfdir , "r") as GFfile:
        trRP = GFfile[hdir + "RP/GF." + dist_form + ".SY.LHT.SAC"][...]
        trTP = GFfile[hdir + "TP/GF." + dist_form + ".SY.LHT.SAC"][...]

    return (trRP, trTP)



def MTrotationE_hdf5(azi, (trRPsy, trTPsy)):
    """
    Rotates the Greens function for the given azimuth (*azi*).

    :param float azi: The station azimuth with respect to the hypocenter.
    :param (trPPsy,  trRRsy, trRTsy,  trTTsy): Output of :py:func:`GFSelectZ`.

    :return: Tuple containin the four greens functions corresponding to the vertical component.
    :rtype: Tuple of six arrays.
    """

    sinp = np.sin(azi)
    cosp = np.cos(azi)
    sin2p = np.sin(2.*azi)
    cos2p = np.cos(2.*azi)

    #Obtaining the rotated synthetics:
    trRRsy_rot =   np.zeros(trRPsy.shape,dtype=np.float32)
    trPPsy_rot =   0.5*sin2p*trTPsy
    trTTsy_rot =  -0.5*sin2p*trTPsy
    trTPsy_rot =   cos2p*trTPsy
    trRTsy_rot =   -sinp*trRPsy
    trRPsy_rot =  cosp*trRPsy

    return [trRRsy_rot, trPPsy_rot,  trTTsy_rot, trTPsy_rot,trRTsy_rot, trRPsy_rot]



def search_best_hdir_hdf5(hypdepth,hdirs):

    depths = [float(dep[1:]) for dep in hdirs]
    BestDirIndex = np.argsort(abs(hypdepth-np.array(depths)))[0]
    hdir_best = hdirs[BestDirIndex] + "/"
    # hdir is the absolute path to the closest deepth.
    return hdir_best



def ltrim(data, starttime, delta):
    '''
    Left trimming similar to obspy but when *data* is a np array.

    :param data: one dimensional array to trim.
    :type data: :py:class:`numpy.array`.
    :param numeric starttime: The amount of time to trim from the front of the trace in seconds.
    :param float delta: Sampling interval in seconds.

    :returns: A view of the relevant subset of *data*.
    '''

    i_of = int(round(starttime/delta))
    if i_of < 0:
        gap = np.empty(abs(i_of))
        gap.fill(data[0])
        data = np.concatenate((gap, data))
        return data
    return data[i_of:]

