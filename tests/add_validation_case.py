#!/usr/bin/env python3

import json
import logging
import os
import sys
from os.path import join, abspath, dirname

from obspy import UTCDateTime

from eatws_skip_client import SKIP, Event
from wphase import runwphase
from validation_cases import result_keys, DATA_DIR, add_case

skip = SKIP('https://skip.eatws.net',
            secret_id='skip-prod-readonly-access')

logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

def prepare_test_case(evid, inventory=None, waveforms=None):
    event: Event = skip.get_event(evid)
    datadir = abspath(DATA_DIR)
    if not event:
        raise ValueError('%s not found in SKIP' % evid)
    eqinfo = dict(
        id=event.id,
        lat=event.latitude,
        lon=event.longitude,
        dep=event.depth_km,
        time=UTCDateTime(event.event_time),
    )
    result = runwphase(
        server='IRIS',
        eqinfo=eqinfo,
        save_waveforms=join(datadir, '%s.mseed' % evid),
        save_inventory=join(datadir, '%s.xml' % evid),
        inventory=inventory,
        waveforms=waveforms,
    )
    MT = result['MomentTensor']
    case = dict(
        id=event.id,
        lat=event.latitude,
        lon=event.longitude,
        dep=event.depth_km,
        time=event.event_time,
        _expected_results={k: MT[k] for k in result_keys},
    )
    print(json.dumps(case, indent=4))
    add_case(case)
    print("This test case has been added to validation_cases.json and test-datasets/.")
    print("To create a new release tarball: "
          "tar czvf ga-wphase-test-datasets.tar.gz test-datasets/ validation_cases.json")

if __name__ == '__main__':
    prepare_test_case(*sys.argv[1:])
