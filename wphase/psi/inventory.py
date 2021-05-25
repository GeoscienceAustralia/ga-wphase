from __future__ import absolute_import, print_function

from future import standard_library
standard_library.install_aliases()
from builtins import range
import pickle as pickle
from collections import Counter, defaultdict
import numpy as np
from obspy.core import UTCDateTime, Stream

try:
    from obspy.clients.fdsn import Client
except ImportError:
    from obspy.fdsn import Client

try:
    from obspy.geodetics import locations2degrees
except ImportError:
    from obspy.core.util.geodetics import locations2degrees

from wphase.psi.taup_fortran import getPtime
import wphase.psi.seismoutils as SU
from wphase.psi.decimate import decimateTo1Hz, CannotDecimate
from wphase import logger

def build_metadata_dict(
        inv,
        target_sampling_rates=[1., 20., 40., 50.],
        target_channels=['BH', 'LH'],
        target_locs=None):
    '''
    Builds a metadata catalog given an ObsPy inventory object.
    The output dictonary is designed to be fully compatible with the
    :py:func:`wprealtime_FFI.wpinv`.

    :param inv: The inventory.
    :type inv: :py:class:`obspy.core.invventory`
    :param target_sampling_rates: The sampling rates to keep metadata for. If *None*, then keep
        all sampling rates.
    :type target_sampling_rates: list or None
    :param target_channels: The channels to keep. If *None*, then keep all channels.
    :type target_channels: list or None
    :param target_locs: The location codes to keep. If *None*, then keep all locations.
    :type target_locs: list or None
    :param outfile: The name of a file to write to catalog to (using :py:func:`cPickle.dump`).
         If *None*, then don't write to file.
    :type outfile: str

    :return: Dictionary containing metadats suitible for passing to
        :py:class:`wphase.psi.wprealtime_FFI.wpinv`.
    :rtype: dict
    '''

    metadata = {}
    failures = []
    def include_channel(cha):
        if target_channels is not None:
            if cha.code[:-1] not in target_channels:
                return False
        if target_sampling_rates is not None:
            if cha.sample_rate not in target_sampling_rates:
                return False
        if target_locs is not None:
            if cha.location_code not in target_locs:
                return False
        return True

    transfer_function_mapping = defaultdict(lambda: "Unknown", {
        "LAPLACE (RADIANS/SECOND)": "A",
        "LAPLACE (HERTZ)": "B",
    })
    for net in inv:
        for sta in net:
            for cha in sta:
                try:
                    if not include_channel(cha):
                        continue
                    trid = '.'.join((net.code, sta.code, cha.location_code, cha.code))
                    resp = cha.response
                    paz = resp.get_paz()

                    metadata[trid] = dict(
                        latitude=cha.latitude,
                        longitude=cha.longitude,
                        elevation=cha.elevation,
                        azimuth=cha.azimuth,
                        dip=cha.dip,
                        sensitivity=resp.instrument_sensitivity.value,
                        sampling_rate=cha.sample_rate,
                        transfer_function=transfer_function_mapping[paz.pz_transfer_function_type],
                        zeros=paz.zeros,
                        poles=paz.poles,
                        gain=paz.normalization_factor,
                    )
                except Exception as e:
                    logger.warning(str(e))
                    failures.append(trid)

    return metadata, failures



def get_waveforms(
        eqinfo,
        META,
        wp_tw_factor=15,
        t_beforeP=1500., # seconds
        t_afterWP=60.,   # seconds
        client=None,
        dist_range=(5., 90.), # degrees
        add_ptime=True,
        bulk_chunk_len=200,
        prune_cutoffs=(1., 2., 3., 4., 5., 5.),
        decimate=True,
        reject_incomplete=False,
        req_times=None,
        waveforms=None):
    '''
    This function will get the waveforms associated with an event to
    perform a WP inversion.

    :param float wp_tw_factor: Defines the Wphase time window. this is used as
        :math:`[t_p, t_p + wp\_tw\_factor \Delta]`.
    :param float dist_range: A pair of floats specifying the min/max epicentral
        distance (in degrees) to be considered; i.e. stations outside of this
        range will be excluded.
    :param bool add_ptime: Should the p arrival time (with respect to the origin
        time) be included in each stations metadata. Including this should
        result in a significant speedup.
    :param int bulk_chunk_len: Maximum number of entries (channels) to be
        included in a single request to a metadata service.
    :param prune_cutoffs: List of pruning (distance) cutoffs passed to
        :py:func:`station_pruningNEZ` :type prune_cutoffs: list of floats
    :param bool decimate: Perform decimation of BH channels. **This should
        always be true when running Wphase**.
    :param dict req_times: Dictionary keyed by channel id containing the start
        and end times of the time window required by Wphase.
    :param client: An Obspy FDSN client.

    :return: If *add_ptime* is *True*, then return a two element tuple containing:
        #. An :py:class:`obspy.core.stream.Stream` containing the data.
        #. A new metadata dictionary containing the p arrival time for each trace.

        If *add_ptime* is *False* return the stream only.
    '''

    hyplat = eqinfo['lat']
    hyplon = eqinfo['lon']
    hypdep = eqinfo['dep']
    otime = eqinfo['time']

    # Obtaining stations within the distance range only. trlist_in_dist_range will contain
    # the ids of channels which are within the specified distance range.
    trlist_in_dist_range = []
    tr_dists_in_range = []
    for trid, stmeta in META.items():
        stlat, stlon = stmeta['latitude'], stmeta['longitude']
        dist = locations2degrees(hyplat, hyplon, stlat, stlon)
        if dist >= dist_range[0] and dist <= dist_range[1]:
            trlist_in_dist_range.append(trid)
            tr_dists_in_range.append(dist)

    # Determining times for the waveform request. req_times will contain time series
    # window for each channel.
    # {
    #   <channel-id>: (<obspy.core.utcdatetime.UTCDateTime>, <obspy.core.utcdatetime.UTCDateTime>)
    # }
    #
    # new_META:
    # {
    #   <channel-id>: original metadata for channel and extra key 'ptime'
    # }
    if req_times is None:
        req_times = {}
        if add_ptime:
            new_META = {}
        for i_trid, trid in enumerate(trlist_in_dist_range):
            stmeta = META[trid]
            dist = tr_dists_in_range[i_trid]
            t_p = getPtime(dist, hypdep)
            t_p_UTC = otime + t_p
            t_wp_end = wp_tw_factor*dist
            t_wp_end_UTC = t_p_UTC + t_wp_end
            t1 = t_p_UTC - t_beforeP
            t2 = t_wp_end_UTC + t_afterWP
            if add_ptime:
                new_META[trid] = stmeta
                new_META[trid]['ptime'] = t_p
            req_times[trid] = [t1, t2]

    # ---------------get waveforms--------------
    # Station pruning
    if prune_cutoffs is not None:
        trlist_in_dist_range = station_pruningNEZ(
            trlist_in_dist_range,
            new_META, prune_cutoffs)

    #reject incomplete traces:
    if reject_incomplete:
        now = UTCDateTime()
        trlist_in_dist_range = [trid for trid in trlist_in_dist_range
            if req_times[trid][1] < now]

    if waveforms:
        # waveforms provided as input, just clean them
        logger.info('%d traces provided as input', len(waveforms))
        st = waveforms
    else:
        # fetch waveforms from server
        logger.info('fetching data from %s', client.base_url)
        st = Stream()

        # Create the subsets for each request
        bulk = []
        for trid in trlist_in_dist_range:
            net, sta, loc, cha = trid.split('.')
            t1, t2 = req_times[trid]
            bulk.append([net, sta, loc, cha, t1, t2])
        bulk_chunks = [bulk[i_chunk:i_chunk + bulk_chunk_len]
                       for i_chunk in range(0, len(bulk), bulk_chunk_len)]

        # make a call for each subset
        # TODO: One might want to do this in parallel.
        for chunk in bulk_chunks:
            # TODO: Do want to try/catch here?
            try:
                st += client.get_waveforms_bulk(chunk)
            except Exception as e:
                logger.error('Problem with request from server %s:\n%s', client.base_url, str(e))
                continue

    # Removing gappy traces (that is channel ids that are repeated)
    st = remove_gappy_traces(st)
    logger.info('%s traces remaining after throwing out gappy ones', len(st))

    # Decimating BH channels. This can be done in parallel.
    if decimate:
        st_B_cha_list = [tr for tr in st if tr.id.split('.')[-1][0] == 'B']
        for tr in st_B_cha_list:
            try:
                decimateTo1Hz(tr)
            except CannotDecimate as e:
                logger.info("Removing trace %s - %s", tr.id, e)
                st.remove(tr)

        # Creating contigous arrays for the traces. This may speed up things later.
        for tr in st:
            tr.data = np.ascontiguousarray(tr.data)

    if add_ptime:
        return st, new_META
    else:
        return st

def remove_gappy_traces(st):
    """Given an obspy Stream, remove any traces with gaps.

    We do this by first merging, then looking for duplicate waveformStreamIDs."""
    try:
        st = st.merge(fill_value='interpolate')
    except Exception as e:
        # If there are traces with the same stream ID but different sampling
        # rates, Stream.merge() throws a raw Exception!
        logger.warning("Failed to merge traces: %s", e)

    trlist_data = [tr.id for tr in st]
    rep_ids = [trid for trid, nrep in list(Counter(trlist_data).items())
               if nrep > 1]
    st = Stream(tr for tr in st if tr.id not in rep_ids)
    return st


def station_pruningNEZ(trlist_pre, META, cutoffs):
    '''
    Station pruning based on distance cutoffs. This version works when more
    one channel per station is present in st (Stream with traces to be pruned).

    This is a wrapper for :py:func:`seismoutils.station_pruning` that works for
    multiple channels per station. It works by creating a list for each **station id**
    rather than each **channel id**.

    :param list trlest_pre: List of channel ids.
    :param dict META: Channel metadata dictionary, keyed by channel id.
    :param cutoffs: List of pruning (distance) cutoffs passed to
        :py:func:`seismoutils.station_pruning`
    :type cutoffs: list of floats

    :return: The the list of remaining channel ids.
    :rtype: list
    '''

    sta_list = []
    trlist_1cha = []
    for trid in trlist_pre:
        sta_code = '.'.join(trid.split('.')[:2])
        if sta_code not in sta_list:
            sta_list.append(sta_code)
            # it does not matter which channel (within the station) we use.
            trlist_1cha.append(trid)

    # We pass the meta data for the channel because it contains the location
    trlist_1cha_short = SU.station_pruning(META, trlist_1cha, cutoffs=cutoffs)
    sta_list_short = [trid.split('.')[:2] for trid in trlist_1cha_short]
    st_short_list = [trid for trid in trlist_pre if trid.split('.')[:2] in sta_list_short]
    return st_short_list
