#!/usr/bin/env python3

import json
import logging
import sys
from os.path import join, abspath
from typing import Optional

from obspy import UTCDateTime

from eatws_skip_client import SKIP, Event
from wphase import runwphase
from tests.validation_cases import result_keys, DATA_DIR, add_case

# get event data from production SKIP
skip = SKIP("https://skip.eatws.net", secret_id="skip-prod-readonly-access")

logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)


def prepare_test_case(evid, inventory=None, waveforms=None):
    event: Optional[Event] = skip.get_event(evid)
    if not event:
        raise ValueError("%s not found in SKIP" % evid)
    eqinfo = dict(
        id=event.id,
        lat=event.latitude,
        lon=event.longitude,
        dep=event.depth_km,
        time=UTCDateTime(event.event_time),
    )
    datadir = abspath(DATA_DIR)
    result = runwphase(
        server="IRIS",
        eqinfo=eqinfo,
        save_waveforms=join(datadir, "%s.mseed" % evid),
        save_inventory=join(datadir, "%s.xml" % evid),
        inventory=inventory,
        waveforms=waveforms,
    )
    MT = result.MomentTensor
    case = dict(
        id=event.id,
        lat=event.latitude,
        lon=event.longitude,
        dep=event.depth_km,
        time=event.event_time,
        _expected_results={k: getattr(MT, k) for k in result_keys},
    )
    if result.QualityParams is not None:
        case["_expected_results"]["azimuthal_gap"] = result.QualityParams.azimuthal_gap
    print(json.dumps(case, indent=4))
    add_case(case)
    print("This test case has been added to validation_cases_v2.json and test-datasets/.")
    print(
        "To create a new release tarball: "
        "tar czvf ga-wphase-test-datasets.tar.gz test-datasets/ validation_cases_v2.json"
    )


if __name__ == "__main__":
    prepare_test_case(*sys.argv[1:])
