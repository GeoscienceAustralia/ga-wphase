from pathlib import Path
import json
import logging
import tarfile
import os
from os.path import join, exists
from urllib.request import urlretrieve
from tempfile import NamedTemporaryFile

import obspy
from obspy.core import UTCDateTime

logger = logging.getLogger("wphase.tests")

DATASETS_URL = 'https://github.com/GeoscienceAustralia/ga-wphase/releases/download/v0.3.2/ga-wphase-test-datasets.tar.gz'
TESTS_DIR = Path(__file__).parent
DOWNLOADED_DIR = TESTS_DIR / "downloaded"
BIG_FILES_DIR = DOWNLOADED_DIR / "test-datasets"

def fetch_datasets():
    """Download and extract the test datasets."""
    tarball_path = TESTS_DIR / "dl.tar.gz"
    logger.warning("Downloading test datasets from %s", DATASETS_URL)
    urlretrieve(DATASETS_URL, tarball_path)
    with tarfile.open(tarball_path, "r:gz") as tarball:
        logger.warning("Extracting test datasets to %s", DOWNLOADED_DIR)
        tarball.extractall(path=DOWNLOADED_DIR)
    os.remove(tarball_path)


def get_dataset(eqinfo):
    """Retrieve a test dataset, either from the local cache directory or from
    the web."""
    evid = eqinfo["id"]
    wfpath = BIG_FILES_DIR / f"{evid}.mseed"
    invpath = BIG_FILES_DIR / f"{evid}.xml"
    if not (wfpath.exists() and invpath.exists()):
        fetch_datasets()
    if not (wfpath.exists() and invpath.exists()):
        raise Exception(f"Dataset {evid} missing even after running fetch_datasets!")
    inventory = obspy.read_inventory(invpath)
    waveforms = obspy.read(wfpath)
    return inventory, waveforms


result_keys = [
    "drlat",
    "drlon",
    "drmag",
    "drdepth",
    "tmtp",
    "tmtt",
    "tmrt",
    "tmrr",
    "tmrp",
    "tmpp",
]

def parse_case(case):
    case['time'] = UTCDateTime(case['time'])
    return case

def dump_case(case):
    case = case.copy()
    if not isinstance(case['time'], str):
        case['time'] = case['time'].strftime("%Y-%m-%dT%H:%M:%SZ")
    return case

def _load_cases():
    with (DOWNLOADED_DIR / "validation_cases.json").open() as fh:
        return [parse_case(x) for x in json.load(fh)]

try:
    cases = _load_cases()
except FileNotFoundError:
    fetch_datasets()
    cases = _load_cases()

def add_case(case):
    cases.append(case)
    with open(join(DOWNLOADED_DIR, "validation_cases.json"), "w") as fh:
        json.dump([dump_case(x) for x in cases], fh, indent=4)
