import numpy.testing
import pytest

from wphase import runwphase
from wphase.psi.model import Event

from validation_cases import cases, get_dataset, result_keys

def assert_allclose(a, b):
    numpy.testing.assert_allclose(a, b, rtol=0.01) # we tolerate a 1% error

@pytest.mark.parametrize("event", cases, ids=lambda c: c['id'])
def test_validity_from_fixed_datasets(event, tmpdir):
    """Ensure that the inversion returns the expected results when tested on
    known input data."""
    # Retrieving data dynamically from IRIS means we can get different inputs
    # and thus different results; so we run a validation test using fixed inputs instead.
    # These datasets are several MB in size, so we don't store them in the repository.
    # They are instead downloaded from this public S3 bucket.
    inventory, waveforms = get_dataset(event)

    eqinfo = Event(
        id=event["id"],
        latitude=event["lat"],
        longitude=event["lon"],
        depth=event["dep"],
        time=event["time"],
    )

    r = runwphase(eqinfo=eqinfo,
                  output_dir=str(tmpdir), output_dir_can_exist=True,
                  waveforms=waveforms, inventory=inventory,
                  make_maps=True, make_plots=True,
                  raise_errors=True)

    MT = r.MomentTensor
    assert MT is not None
    assert MT.drmagt == 'Mww'

    EXP = event['_expected_results']
    assert_allclose([EXP[k] for k in result_keys],
                    [getattr(MT, k) for k in result_keys])
    ag = EXP.get("azimuthal_gap")
    if ag is not None:
        assert r.QualityParams is not None
        assert_allclose(r.QualityParams.azimuthal_gap, ag)
