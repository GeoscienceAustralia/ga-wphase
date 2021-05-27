import logging
import tarfile
import os
from os.path import join, exists
from urllib.request import urlretrieve
from tempfile import NamedTemporaryFile

import obspy
from obspy.core import UTCDateTime

logger = logging.getLogger("wphase.tests")

DATASETS_URL = 'https://github.com/GeoscienceAustralia/ga-wphase/releases/download/v0.1/ga-wphase-test-datasets-v0.1.tar.gz'
TARBALL_SUBDIR = 'test-datasets'
TESTS_DIR = os.path.dirname(__file__)

def fetch_datasets(url):
    """Download and extract the test datasets."""
    tarball_path = join(TESTS_DIR, 'dl.tar.gz')
    logger.warning("Downloading test datasets from %s", DATASETS_URL)
    urlretrieve(url, tarball_path)
    with tarfile.open(tarball_path, "r:gz") as tarball:
        logger.warning("Extracting test datasets to %s", TESTS_DIR)
        tarball.extractall(path=TESTS_DIR)
    os.remove(tarball_path)


def get_dataset(eqinfo):
    """Retrieve a test dataset, either from the local cache directory or from
    the web."""
    evid = eqinfo["id"]
    wfpath = join(TESTS_DIR, TARBALL_SUBDIR, "{}.mseed".format(evid))
    invpath = join(TESTS_DIR, TARBALL_SUBDIR, "{}.xml".format(evid))
    if not (exists(wfpath) and exists(invpath)):
        fetch_datasets(DATASETS_URL)
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

cases = [
    {
        "id": "ga2016rhegtj",
        "lat": -37.228,
        "lon": 178.755,
        "dep": 40.5,
        "time": UTCDateTime("2016-09-01T16:38:00Z"),
        "_expected_results": {
            "drlat": -37.028000000000013,
            "drlon": 178.95500000000004,
            "drmag": 7.0738536889439168,
            "drdepth": 30.5,
            "tmtp": 1.1570247877704712e19,
            "tmtt": -6.652347385786327e18,
            "tmrt": 2.2844251194450883e19,
            "tmrr": -1.9274432292250579e19,
            "tmrp": 3.7947744865955242e19,
            "tmpp": 2.5926779678036906e19,
        },
    },
    {
        "id": "ga2019kgtpbq",
        "lat": -5.806561947,
        "lon": -75.21974182,
        "time": UTCDateTime("2019-05-26T07:41:14Z"),
        "dep": 115.8131256,
        "_expected_results": {
            "drlat": -5.6065619469999959,
            "drlon": -75.019741820000007,
            "drmag": 8.0595872064415115,
            "drdepth": 120.5,
            "tmtp": -3.007618e+20,
            "tmtt": -1.061544e+20,
            "tmrt": 2.955757e+20,
            "tmrr": -1.262713e+21,
            "tmrp": -6.980889e+20,
            "tmpp": 1.368867e+21,
        },
    },
    {
        "id": "ga2020nzrmhv",
        "lat": -7.849661827,
        "lon": 147.7071838,
        "time": UTCDateTime("2020-07-17T02:50:22Z"),
        "dep": 81.35980988,
        "_expected_results": {
            "drlat": -7.6496618269999974,
            "drlon": 147.90718380000004,
            "drmag": 7.031759280578866,
            "drdepth": 100.5,
            "tmtp": 7.3971395108235151e18,
            "tmtt": -2.2987138457087701e19,
            "tmrt": -1.482707570357497e19,
            "tmrr": -6.7682667799839885e17,
            "tmrp": 3.398088023337746e19,
            "tmpp": 2.36639651350861e19,
        },
    },
    {
        "id": "ga2020twkduz",
        "lat": -6.114739418,
        "lon": 146.1887207,
        "time": UTCDateTime("2020-10-08T07:35:32Z"),
        "dep": 107.0434341,
        "_expected_results": {
            "drlat": -5.9147394179999964,
            "drlon": 146.38872070000005,
            "drmag": 6.4892603760544638,
            "drdepth": 90.5,
            "tmtp": -2.9266909677544668e18,
            "tmtt": 4.7587623734491853e18,
            "tmrt": -3.3093074468640128e18,
            "tmrr": 7.461857380934272e17,
            "tmrp": 5.1437552271745523e17,
            "tmpp": -5.504948111542612e18,
        },
    },
    {
        "id": "ga2019mhqctz",
        "lat": -6.424966812,
        "lon": 129.1710968,
        "time": UTCDateTime("2019-06-24T02:53:39Z"),
        "dep": 210.3308105,
        "_expected_results": {
            "drlat": -6.2249668119999981,
            "drlon": 129.37109679999998,
            "drmag": 7.2421266535481816,
            "drdepth": 200.5,
            "tmtp": 5.927167106444483e19,
            "tmtt": -7.1178408071732183e19,
            "tmrt": -3.1365589411401138e19,
            "tmrr": 2.7839662569089118e19,
            "tmrp": 9.1634140771453071e18,
            "tmpp": 4.3338745502643061e19,
        },
    },
    {
        "id": "ga2011excpfy",
        "lat": 38.27,
        "lon": 142.51,
        "time": UTCDateTime("2011-03-11T05:46:21Z"),
        "dep": 10,
        "_expected_results": {
            "drlat": 37.67,
            "drlon": 142.71,
            "drmag": 9.096786,
            "drdepth": 13.5,
            "tmtp": -6.476268e21,
            "tmtt": -1.661770e21,
            "tmrt": 1.248122e22,
            "tmrr": 1.719292e22,
            "tmrp": 5.123794e22,
            "tmpp": -1.553115e22,
        },
    },
    {
        "id": "ga2021cwlolj",
        "lat": -22.77105522,
        "lon": 171.786438,
        "time": UTCDateTime("2021-02-10T13:19:58Z"),
        "dep": 22.85298347,
        "_expected_results": {
            "drlat": -22.571055219999998,
            "drlon": 171.18643800000004,
            "drmag": 7.839797036936293,
            "drdepth": 15.5,
            "tmtp": 3.5182321594404639e19,
            "tmtt": -2.8019089697167357e20,
            "tmrt": 6.446917682329134e20,
            "tmrr": 3.0545559224700869e20,
            "tmrp": 1.4485417735880932e20,
            "tmpp": -2.5264695275335123e19,
        },
    },
]
