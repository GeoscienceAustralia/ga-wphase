#!/usr/bin/env python3
"""
Re-run wphase against existing validation datasets and update
_expected_results in validation_cases.json.

This script:
- Uses the locally stored test datasets (mseed + inventory XML)
- Does NOT fetch data from IRIS / SKIP
- Recomputes expected results for all cases
"""
from tempfile import TemporaryDirectory

import json
import logging
from copy import deepcopy
from os.path import join
from pathlib import Path

from wphase import runwphase
from wphase.psi.model import Event

from tests.validation_cases import (
    cases,
    get_dataset,
    result_keys,
    dump_case,
    DOWNLOADED_DIR,
)

logging.basicConfig(level=logging.INFO)
log = logging.getLogger("wphase.update_validation")


def recompute_case(case):
    log.info("Recomputing expected results for %s", case["id"])

    inventory, waveforms = get_dataset(case)

    eqinfo = Event(
        id=case["id"],
        latitude=case["lat"],
        longitude=case["lon"],
        depth=case["dep"],
        time=case["time"],
    )

    with TemporaryDirectory() as tdir:
        result = runwphase(
            eqinfo=eqinfo,
            waveforms=waveforms,
            inventory=inventory,
            raise_errors=True,
            make_maps=False,
            make_plots=False,
            output_dir=tdir,
            output_dir_can_exist=True,
        )

    MT = result.MomentTensor
    if MT is None:
        raise RuntimeError(f"No MomentTensor produced for {case['id']}")

    expected = {k: getattr(MT, k) for k in result_keys}

    if result.QualityParams is not None:
        expected["azimuthal_gap"] = result.QualityParams.azimuthal_gap

    return {**case, "_expected_results": expected}


def main():
    updated_cases = [recompute_case(c) for c in cases]

    outpath = Path(DOWNLOADED_DIR) / "validation_cases.json"
    log.info("Writing updated expected results to %s", outpath)

    with outpath.open("w") as fh:
        json.dump([dump_case(c) for c in updated_cases], fh, indent=4)

    log.info("Done. validation_cases.json updated.")


if __name__ == "__main__":
    main()
