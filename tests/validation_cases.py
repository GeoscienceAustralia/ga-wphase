import json
import logging
import tarfile
import os
from os.path import join, exists, dirname
from urllib.request import urlretrieve
from tempfile import NamedTemporaryFile

import obspy
from obspy.core import UTCDateTime

logger = logging.getLogger("wphase.tests")

DATASETS_URL = 'https://github.com/GeoscienceAustralia/ga-wphase/releases/download/v0.3/ga-wphase-test-datasets.tar.gz'
DATA_DIR = 'test-datasets'
TESTS_DIR = dirname(__file__)

def fetch_datasets():
    """Download and extract the test datasets."""
    tarball_path = join(TESTS_DIR, 'dl.tar.gz')
    logger.warning("Downloading test datasets from %s", DATASETS_URL)
    urlretrieve(DATASETS_URL, tarball_path)
    with tarfile.open(tarball_path, "r:gz") as tarball:
        logger.warning("Extracting test datasets to %s", TESTS_DIR)
        tarball.extractall(path=TESTS_DIR)
    os.remove(tarball_path)


def get_dataset(eqinfo):
    """Retrieve a test dataset, either from the local cache directory or from
    the web."""
    evid = eqinfo["id"]
    wfpath = join(TESTS_DIR, DATA_DIR, "{}.mseed".format(evid))
    invpath = join(TESTS_DIR, DATA_DIR, "{}.xml".format(evid))
    if not (exists(wfpath) and exists(invpath)):
        fetch_datasets()
    if not (exists(wfpath) and exists(invpath)):
        raise Exception("Dataset {evid} missing even after running fetch_datasets!")
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
    with open(join(TESTS_DIR, "validation_cases.json")) as fh:
        return [parse_case(x) for x in json.load(fh)]

try:
    cases = _load_cases()
except FileNotFoundError:
    fetch_datasets()
    cases = _load_cases()

def add_case(case):
    cases.append(case)
    with open(join(TESTS_DIR, "validation_cases.json"), "w") as fh:
        json.dump([dump_case(x) for x in cases], fh, indent=4)
