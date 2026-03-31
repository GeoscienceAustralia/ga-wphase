import numpy.testing
import pytest

from wphase.psi.seismoutils import get_azimuths
from wphase.psi.core import azimuthal_gap

GAP_CASES = [
    # Check degenerate cases:
    ([], 360),
    ([20], 360),
    # Check wrap-around:
    ([20, 30], 350),
    ([20, 230], 210),
    # Typical cases:
    ([20, 50, 70, 90, 115, 140, 165, 190, 215, 230, 255, 280, 300, 320, 340, 360], 30),
    ([70, 90, 110, 130, 150, 190, 210, 240, 245, 260, 280, 300, 320, 340], 90),
    # Check order doesn't matter:
    ([90, 110, 130, 150, 190, 20, 210, 240, 245, 260, 70, 280, 300, 320, 340, 360], 50),
]

@pytest.mark.parametrize("case", GAP_CASES)
def test_azimuthal_gap(case):
    azis, gap = case
    assert azimuthal_gap(azis) == gap

def test_get_azimuths():
    origin = (70, 140)
    metadata = {
        "E": {"latitude": 70, "longitude": 150},
        "W": {"latitude": 70, "longitude": 130},
        "N": {"latitude": 80, "longitude": 140},
        "N2": {"latitude": 70, "longitude": -40},
        "S": {"latitude": 0, "longitude": 140},
    }
    keys = ["E", "W", "N", "N2", "S"]
    azis = get_azimuths(metadata, keys, origin)
    numpy.testing.assert_allclose(azis, [
       85.30014199486736,
       274.6998580051326,
       0.0,
       360.0,
       180.0,
    ], rtol=1e-9, atol=1e-7)
