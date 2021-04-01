from future import standard_library

standard_library.install_aliases()

import os
from builtins import str

import numpy.testing
import pytest
from wphase import runwphase

from validation_cases import cases, get_dataset, result_keys


def assert_allclose(a, b):
    numpy.testing.assert_allclose(a, b, rtol=1e-3)

@pytest.mark.parametrize("eqinfo", cases)
def test_validity_from_fixed_datasets(eqinfo, tmpdir):
    """Ensure that the inversion returns the expected results when tested on
    known input data."""
    # Retrieving data dynamically from IRIS means we can get different inputs
    # and thus different results; so we run a validation test using fixed inputs instead.
    # These datasets are several MB in size, so we don't store them in the repository.
    # They are instead downloaded from this public S3 bucket.
    inventory, waveforms = get_dataset(eqinfo)

    r = runwphase(eqinfo=eqinfo,
                  output_dir=str(tmpdir), output_dir_can_exist=True,
                  waveforms=waveforms, inventory=inventory,
                  make_maps=False, make_plots=False,
                  raise_errors=True)

    assert 'MomentTensor' in r

    MT = r['MomentTensor']
    assert MT['drmagt'] == 'Mww'

    EXP = eqinfo['_expected_results']
    assert_allclose([EXP[k] for k in result_keys],
                    [MT[k] for k in result_keys])
