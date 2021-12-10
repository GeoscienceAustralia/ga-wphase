"""Implements runwphase, which acquires data from FDSN and runs the W-Phase
inversion, including post-processing."""
import logging
import os
import pickle
from functools import partial
from traceback import format_exc

import obspy
from obspy.clients.fdsn import Client
from obspy.core.inventory import Inventory
from obspy.core.utcdatetime import UTCDateTime

from wphase import settings
from wphase.data_acquisition import build_metadata_dict, get_waveforms
from wphase.psi.core import wpinv
from wphase.psi.exceptions import InversionError
from wphase.psi.model import Event, WPhaseResult
from wphase.wputils import (NoProfiler, OutputDict, WPInvProfiler,
                            post_process_wpinv)

logger = logging.getLogger(__name__)


class ArrayLogger(logging.Handler):
    """Handler that captures log mesages and stores them in an array."""
    def __init__(self, *args, **kwargs):
        super(ArrayLogger, self).__init__(*args, **kwargs)
        self.messages = []

    def emit(self, record):
        self.messages.append(self.format(record).strip())


class LogCapture(object):
    """Context manager to capture logs."""
    def __init__(self, logger, level):
        self.logger = logger
        self.level = level
        self.handler = ArrayLogger(level=level)

    def __enter__(self):
        self.logger.addHandler(self.handler)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.logger.removeHandler(self.handler)
        self.messages = self.handler.messages


def load_metadata(
        client: Client,
        eqinfo: Event,
        dist_range,
        networks,
        inventory=None,
        save_path=None,
        t_before_origin=3600.,
        t_after_origin=3600.):

    """
    :param float t_after_origin: The end of the time window around
        *eqinfo['time']* (in seconds) in which the station must exist. This
        should not be confused with the parameter *t_afterWP* used elsewhere,
        which is used in a more complicated calculation which requires the
        location of the station, which we get in this function. This is only
        used to filter the request for the inventory and hence only need be
        very rough (within the week would probably be equally sensible).
    """
    if inventory:
        if isinstance(inventory, str):
            inventory = obspy.read_inventory(inventory)
        logger.info('Loaded provided inventory')
        inv = inventory
    else:
        def caller_maker(depth=0, **kwargs):
            if 'network' in kwargs and kwargs['network'].upper() == 'ALL':
                kwargs.pop('network')

            base_call = {
                "level"    : 'response',
                "channel"  : 'BH?',
                "latitude" : eqinfo.latitude,
                "longitude": eqinfo.longitude,
                "minradius": dist_range[0],
                "maxradius": dist_range[1],
                "starttime": eqinfo.time - t_before_origin,
                "endtime"  : eqinfo.time + t_after_origin}

            base_call.update(kwargs)

            def make_call(**kwargs):
                args = base_call.copy()
                args.update(kwargs)
                return client.get_stations(**args)

            logger.info('Retrieving metadata from server %s', client.base_url)
            args = dict()
            if depth >= 1:
                args["network"] = net
            if depth >= 2:
                args["station"] = sta
            if depth >= 3:
                args["channel"] = cha

            return partial(make_call, **args)

        try:
            # first, try and get everything
            inv = caller_maker(network=networks)()

        except Exception:
            # ... that didn't work
            nets = caller_maker(network=networks, level='network')()
            inv = Inventory([], None)

            # try by network
            call1 = caller_maker(1)
            for net in nets:
                try:
                    inv += call1(net.code)
                except Exception:
                    # ... by station
                    stas = caller_maker(network=net.code, level='station')()
                    call2 = caller_maker(2)
                    for sta in stas[0]:
                        try:
                            inv += call2(net.code, sta.code)
                        except Exception:
                            # ... by channel
                            chans = caller_maker(network=net.code, station=sta.code, level='channel')()
                            call3 = caller_maker(3)
                            for chan in chans[0][0]:
                                try:
                                    inv += call3(net.code, sta.code, chan.code)
                                except Exception:
                                    # ... skip the channel
                                    # TODO: log that this has happenned
                                    pass

    if save_path:
        logger.info("Saving inventory in %s", save_path)
        inv.write(save_path, format='STATIONXML')
    return build_metadata_dict(inv)



def runwphase(
        output_dir,
        server = None,
        greens_functions_dir = settings.GREENS_FUNCTIONS,
        n_workers_in_pool = settings.WORKER_COUNT,
        processing_level = 3,
        networks = 'II,IU',
        eqinfo: Event = None,
        wp_tw_factor = 15,
        t_beforeP = 1500.,
        t_afterWP = 60.,
        dist_range = [5.,90.],
        add_ptime = True,
        bulk_chunk_len = 200,
        prune_cutoffs = None,
        use_only_z_components = True,
        user=None,
        password=None,
        inventory=None,
        waveforms=None,
        pickle_inputs=False,
        make_maps=True,
        make_plots=True,
        raise_errors=False,
        save_waveforms=None,
        save_inventory=None,
        **kwargs) -> WPhaseResult:

    """
    Run wphase.

    :param output_dir: Full file path to the output directory. **DO NOT USE
        RELATIVE PATHS**.
    :param greens_functions_dir: The Green data Directory.
    :param n_workers_in_pool: Number of processors to use, (default
        :py:data:`wphase.settings.WORKER_COUNT`) specifies as many as is
        reasonable'.
    :param processing_level: Processing level.
    """

    if server is None and (inventory is None or waveforms is None):
        raise ValueError("If not providing server, you must provide inventory and waveforms.")

    meta_t_p = {}

    if eqinfo is None:
        raise ValueError('eqinfo cannot be None')

    if server is not None:
        client = Client(server, user=user, password=password)
    else:
        client = None

    # get the metadata for the event
    metadata, failures = load_metadata(
            client,
            eqinfo,
            dist_range,
            networks,
            inventory=inventory,
            save_path=save_inventory)

    if output_dir and failures:
        with open(os.path.join(output_dir, 'inv.errs'), 'w') as err_out:
            err_out.write('\n'.join(failures))

    if not metadata:
        raise Exception('no metadata available for: {}'.format(eqinfo))

    if output_dir and pickle_inputs:
        with open(os.path.join(output_dir, 'inv.pkl'), 'wb') as inv_out:
            pickle.dump(metadata, inv_out)

    # This is replaced by the return value of wpinv unless there's a fatal error:
    wphase_output = WPhaseResult(Event=eqinfo)

    try:
        # load the data for from the appropriate server
        streams, meta_t_p_ = get_waveforms(
            eqinfo,
            metadata,
            wp_tw_factor = wp_tw_factor,
            t_beforeP = t_beforeP,
            t_afterWP = t_afterWP,
            client = client,
            dist_range = dist_range,
            add_ptime = add_ptime,
            bulk_chunk_len = bulk_chunk_len,
            prune_cutoffs = prune_cutoffs,
            waveforms = waveforms,
            save_path = save_waveforms,
        )

        if use_only_z_components:
            streams = streams.select(component = 'Z')

            logger.info('%d traces remaining after restricting to Z', len(streams))

        meta_t_p.update(meta_t_p_)

        if output_dir and pickle_inputs:
            streams_pickle_file = os.path.join(output_dir, 'streams.pkl')
            with open(streams_pickle_file, 'wb') as pkle:
                pickle.dump((meta_t_p, streams), pkle)

        # do and post-process the inversion
        profiler = WPInvProfiler(output_dir) if settings.PROFILE else NoProfiler()
        with profiler:
            with LogCapture(logger, logging.WARNING) as capture:
                try:
                    wphase_output = wpinv(
                        streams,
                        meta_t_p,
                        eqinfo,
                        greens_functions_dir,
                        processes = n_workers_in_pool,
                        OL = processing_level)
                except InversionError as e:
                    wphase_output.add_warning(e)

            for message in capture.messages:
                wphase_output.add_warning(message)

            try:
                post_process_wpinv(
                    output = wphase_output,
                    WPOL = processing_level,
                    working_dir = output_dir,
                    eqinfo = eqinfo,
                    metadata = meta_t_p,
                    make_maps=output_dir and make_maps,
                    make_plots=output_dir and make_plots)
            except Exception as e:
                if raise_errors:
                    raise
                wphase_output.add_warning("Error during post-processing. %s" % format_exc())
        wphase_output.WPInvProfile = profiler.html

    except Exception as e:
        if raise_errors:
            raise
        wphase_output.Error = str(e)
        wphase_output.StackTrace = format_exc()

    return wphase_output
