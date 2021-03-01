import pytest
from obspy import UTCDateTime

from wphase import runwphase

ga2016rhegtj = {
    "id": "ga2016rhegtj",
    "lat": -37.228,
    "lon": 178.755,
    "dep": 40.5,
    "time": UTCDateTime("2016-09-01T16:38:00Z"),
}

@pytest.mark.slow
def test_ga2016rhegtj_with_iris_fdsn(tmpdir):
    """Test that runwphase can succeed when sourcing data from a public FDSN
    server."""
    # We run one test with IRIS, but don't validate the exact results.
    r = runwphase(eqinfo=ga2016rhegtj, server='IRIS',
                  output_dir=str(tmpdir), output_dir_can_exist=True,
                  make_maps=False)
    assert isinstance(r['MomentTensor']['drmag'], float)
