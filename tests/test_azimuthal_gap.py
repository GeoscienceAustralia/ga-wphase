import pytest
from wphase.psi.core import azimuthal_gap

CASES = [
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

@pytest.mark.parametrize("case", CASES)
def test_azimuthal_gap(case):
    azis, gap = case
    assert azimuthal_gap(azis) == gap
