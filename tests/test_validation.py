import numpy.testing as nptest
from pandas.testing import assert_series_equal
import pandas as pd
import pytest

from wphase import runwphase

from validation_cases import cases, ValidationCase, mt_as_series


def assert_allclose(a, b):
    nptest.assert_allclose(a, b, rtol=0.01)  # we tolerate a 1% error


@pytest.mark.parametrize("case", cases, ids=lambda c: c.Event.id)
def test_validity_from_fixed_datasets(case: ValidationCase, tmpdir):
    """Ensure that the inversion returns the expected results when tested on
    known input data."""
    # Retrieving data dynamically from IRIS means we can get different inputs
    # and thus different results; so we run a validation test using fixed inputs instead.
    # These datasets are several hundred MB in size, so we don't store them in the repository.
    # They are instead fetched from a file attached to a github release - see
    # validation_cases.py for details.
    inventory, waveforms = case.get_dataset()

    r = runwphase(
        eqinfo=case.Event,
        output_dir=str(tmpdir),
        output_dir_can_exist=True,
        waveforms=waveforms,
        inventory=inventory,
        make_maps=True,
        make_plots=True,
        raise_errors=True,
    )

    assert r.MomentTensor is not None
    assert r.MomentTensor.drmagt == "Mww"

    assert_series_equal(
        mt_as_series(r.MomentTensor),
        mt_as_series(case.MomentTensor),
        rtol=1e-2,
    )
    ag = case.azimuthal_gap
    if ag is not None:
        assert r.QualityParams is not None
        assert_allclose(r.QualityParams.azimuthal_gap, ag)
