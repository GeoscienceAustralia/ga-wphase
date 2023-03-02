import logging
import tarfile
from pathlib import Path
from typing import Optional
from urllib.request import urlretrieve

import obspy
import pandas as pd

from wphase.psi.model import AntelopeMomentTensor, Data, Event, WPhaseResult

logger = logging.getLogger("wphase.tests")

DATASETS_URL = 'https://github.com/GeoscienceAustralia/ga-wphase/releases/download/v0.3.1/ga-wphase-test-datasets.tar.gz'
TESTS_DIR = Path(__file__).parent
DATA_DIR = TESTS_DIR / "test-datasets"


class TestEvent(Event):
    id: str # to make this not null


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
class ExpectedMT(AntelopeMomentTensor):
    """The expected results we check don't include all of the required fields
    in AntelopeMomentTensor, so we have this subclass to allow loading the
    partial records."""
    scm: Optional[float] = None
    drmagt: Optional[float] = None
    str1: Optional[float] = None
    dip1: Optional[float] = None
    rake1: Optional[float] = None
    str2: Optional[float] = None
    dip2: Optional[float] = None
    rake2: Optional[float] = None
    auth: Optional[str] = None

    @staticmethod
    def from_result(result: WPhaseResult):
        mt = result.MomentTensor
        assert mt is not None
        return ExpectedMT(
            drlat=mt.drlat,
            drlon=mt.drlon,
            drmag=mt.drmag,
            drdepth=mt.drdepth,
            tmtp=mt.tmtp,
            tmtt=mt.tmtt,
            tmrt=mt.tmrt,
            tmrr=mt.tmrr,
            tmrp=mt.tmrp,
            tmpp=mt.tmpp,
        )


def mt_as_series(mt: AntelopeMomentTensor):
    return pd.Series({k: getattr(mt, k) for k in result_keys})


class ValidationCase(Data):
    Event: TestEvent
    MomentTensor: ExpectedMT
    azimuthal_gap: Optional[float] = None

    def get_dataset(self):
        """Retrieve a test dataset, either from the local cache directory or from
        the web."""
        evid = self.Event.id
        wfpath = DATA_DIR / f"{evid}.mseed"
        invpath = DATA_DIR / f"{evid}.xml"
        if not (wfpath.is_file() and invpath.is_file()):
            fetch_datasets()
        if not (wfpath.is_file() and invpath.is_file()):
            raise Exception(f"Dataset {evid} missing even after running fetch_datasets!")
        inventory = obspy.read_inventory(invpath)
        waveforms = obspy.read(wfpath)
        return inventory, waveforms

    @staticmethod
    def from_result(result: WPhaseResult):
        e = result.Event
        return ValidationCase(
            Event=TestEvent(**dict(e)),
            MomentTensor=ExpectedMT.from_result(result),
            azimuthal_gap=result.QualityParams and result.QualityParams.azimuthal_gap or None,
        )

def fetch_datasets():
    """Download and extract the test datasets."""
    tarball_path = TESTS_DIR / "dl.tar.gz"
    logger.warning("Downloading test datasets from %s", DATASETS_URL)
    urlretrieve(DATASETS_URL, tarball_path)
    with tarfile.open(tarball_path, "r:gz") as tarball:
        logger.warning("Extracting test datasets to %s", TESTS_DIR)
        tarball.extractall(path=TESTS_DIR)
    tarball_path.unlink()


class CasesNotFound(Exception):
    pass


def _load_cases():
    jfiles = [file for file in DATA_DIR.iterdir() if file.suffix == ".json"]
    if not jfiles:
        raise CasesNotFound()
    return [ValidationCase.parse_file(file) for file in jfiles]


try:
    cases = _load_cases()
except (FileNotFoundError, CasesNotFound):
    fetch_datasets()
    cases = _load_cases()


def add_case(case: ValidationCase):
    cases.append(case)
    (DATA_DIR / case.id).write_text(case.json(indent=4))
