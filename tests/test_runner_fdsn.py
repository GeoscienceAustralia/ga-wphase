import numpy as np
from obspy.core import UTCDateTime

from wphase._runner_fdsn import runwphase

def close(a, b):
    return np.isclose(a, b, rtol=1e-3)

def test_ga2016rhegtj_with_iris_fdsn(tmpdir):
    # A mag ~7 near NZ
    ga2016rhegtj = {'lat': -37.228,
                    'lon': 178.755,
                    'dep': 40.5,
                    'time': UTCDateTime("2016-09-01T16:38:00Z")}

    r = runwphase(eqinfo=ga2016rhegtj, server='IRIS', output_dir=str(tmpdir))
    MT = r['MomentTensor']
    assert MT['drmagt'] == 'Mww'
    assert close(MT['drlon'], 178.955)
    assert close(MT['drlat'], -37.028)
    assert close(MT['drmag'], 7.075)
    assert close(MT['str1'], 353.57)
    assert close(MT['str2'], 210.60)
    assert close(MT['dip1'], 16.964)
    assert close(MT['dip2'], 76.312)
    assert close(MT['rake1'], 234.19)
    assert close(MT['rake2'], -79.88)
