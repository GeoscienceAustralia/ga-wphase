import logging
import sys
import json
from os.path import join, abspath

import obspy
from obspy import UTCDateTime

from tests.validation_cases import cases, add_case, result_keys, TESTS_DIR
from wphase import runwphase
from wphase.psi.model import Event

logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

def update_case(event):
    eqinfo = Event(
        id=event["id"],
        latitude=event["lat"],
        longitude=event["lon"],
        depth=event["dep"],
        time=event["time"],
    )

    wfpath = join(TESTS_DIR, 'test-datasets-v2', "{}.mseed".format(eqinfo.id))
    invpath = join(TESTS_DIR, 'test-datasets-v2', "{}.xml".format(eqinfo.id))
    inventory = obspy.read_inventory(invpath)
    waveforms = obspy.read(wfpath)

    result = runwphase(
        server="IRIS",
        eqinfo=eqinfo,
        inventory=inventory,
        waveforms=waveforms,
        raise_errors=True,
    )
    print(result)
    MT = result.MomentTensor
    case = {
        **event,
        "_expected_results": {k: getattr(MT, k) for k in result_keys},
    }
    if not isinstance(case['time'], str):
        case['time'] = case['time'].strftime("%Y-%m-%dT%H:%M:%SZ")

    if result.QualityParams is not None:
        case["_expected_results"]["azimuthal_gap"] = result.QualityParams.azimuthal_gap
    print(json.dumps(case, indent=4))
    add_case(case)


if __name__ == "__main__":
    for event in cases:
        update_case(event)
    print("Test cases have been updated in validation_cases_v2.json.")
    print(
        "To create a new release tarball: "
        "tar czvf ga-wphase-test-datasets.tar.gz test-datasets/ validation_cases_v2.json"
    )
