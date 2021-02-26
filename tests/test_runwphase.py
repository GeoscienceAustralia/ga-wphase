import os
try:
    from urllib.request import urlretrieve
except ImportError:
    from urllib import urlretrieve

import numpy.testing
from obspy.core import UTCDateTime
import obspy

from wphase import runwphase

def assert_allclose(a, b):
    numpy.testing.assert_allclose(a, b, rtol=1e-3)

DIR = os.path.join(os.path.dirname(__file__), 'temp')

# A mag ~7 near NZ
ga2016rhegtj = {'lat': -37.228,
                'lon': 178.755,
                'dep': 40.5,
                'time': UTCDateTime("2016-09-01T16:38:00Z")}

def test_ga2016rhegtj_from_fixed_datasets(tmpdir):
    # Retrieving data dynamically from IRIS means we can get different inputs
    # and thus different results; so we run a validation test using fixed inputs instead.
    # These datasets are several MB in size, so we don't store them in the repository.
    # They are instead downloaded from this public S3 bucket.
    SERVER = 'http://eatws-public.s3.amazonaws.com'
    wfpath = '/ga2016rhegtj.mseed'
    invpath = '/ga2016rhegtj.xml'
    if not os.path.isdir(DIR):
        os.mkdir(DIR)
    if not os.path.exists(DIR + wfpath):
        urlretrieve(SERVER + wfpath, DIR + wfpath)
    if not os.path.exists(DIR + invpath):
        urlretrieve(SERVER + invpath, DIR + invpath)

    waveforms = obspy.read(DIR + wfpath)
    inventory = obspy.read_inventory(DIR + invpath)

    r = runwphase(eqinfo=ga2016rhegtj,
                  output_dir=str(tmpdir), output_dir_can_exist=True,
                  waveforms=waveforms, inventory=inventory,
                  make_maps=False)

    assert 'MomentTensor' in r

    MT = r['MomentTensor']
    assert MT['drmagt'] == 'Mww'

    expected = (
        ('drlat', -37.028),
        ('drlon', 178.955),
        ('drmag', 7.07610),
        ('drdepth', 30.5),
        ('tmtp', 1.157025e+19),
        ('tmtt', -6.652347e+18),
        ('tmrt', 2.284425e+19),
        ('tmrr', -1.927443e+19),
        ('tmrp', 3.794774e+19),
        ('tmpp', 2.592678e+19),
    )
    assert_allclose([v for k, v in expected], [MT[k] for k, v in expected])

def test_ga2016rhegtj_with_iris_fdsn(tmpdir):
    # We run the test with IRIS, but don't validate the exact results.
    r = runwphase(eqinfo=ga2016rhegtj, server='IRIS',
                  output_dir=str(tmpdir), output_dir_can_exist=True,
                  make_maps=False)
    assert isinstance(r['MomentTensor']['drmag'], float)
