#!/usr/bin/env python3

import logging
import sys

from obspy import UTCDateTime

from eatws_skip_client import SKIP
from wphase import runwphase
from wphase.psi.model import Event
from validation_cases import ExpectedMT, ValidationCase, result_keys, DATA_DIR, add_case

# get event data from production SKIP
skip = SKIP('https://skip.eatws.net',
            secret_id='skip-prod-readonly-access')

logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

def prepare_test_case(evid, inventory=None, waveforms=None):
    event = skip.get_event(evid)
    if not event:
        raise ValueError('%s not found in SKIP' % evid)
    eqinfo = Event(
        id=event.id,
        latitude=event.latitude,
        longitude=event.longitude,
        depth=event.depth_km,
        time=UTCDateTime(event.event_time),
    )
    result = runwphase(
        server='IRIS',
        eqinfo=eqinfo,
        save_waveforms=DATA_DIR / f"{evid}.mseed",
        save_inventory=DATA_DIR / f"{evid}.xml",
        inventory=inventory,
        waveforms=waveforms,
        make_maps=False,
        make_plots=False,
    )
    case = ValidationCase.from_result(result)
    print(case.json(indent=4))
    add_case(case)
    print("This test case has been added to test-datasets/.")
    print("To create a new release tarball: "
          "tar czvf ga-wphase-test-datasets.tar.gz test-datasets/")

if __name__ == '__main__':
    prepare_test_case(*sys.argv[1:])
