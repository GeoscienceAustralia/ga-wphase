from functools import partial
from wphase.psi.model import WPhaseResult
from wphase.seiscomp import createObjects, writeSCML
from seiscomp.core import Time
from seiscomp.datamodel import (
    AUTOMATIC,
    CENTROID,
    FocalMechanism,
    MomentTensor,
    NodalPlane,
    NodalPlanes,
    Origin,
    Magnitude,
    Tensor,
)
import numpy.testing

assert_almost_equal = partial(numpy.testing.assert_allclose, rtol=1e-3)


def test_create_objects(tests_dir):
    sample = WPhaseResult.parse_file(tests_dir / "ga2023fsbydd-result.json")
    objects = createObjects(sample, "GA")
    fm = objects["focalMechanism"]
    do = objects["derivedOrigin"]
    mm = objects["momentMagnitude"]
    assert isinstance(fm, FocalMechanism)
    assert isinstance(do, Origin)
    assert isinstance(mm, Magnitude)

    assert fm.creationInfo().author() == "GA W-phase"
    assert fm.creationInfo().agencyID() == "GA"
    assert fm.methodID() == "wphase"
    assert fm.evaluationMode() == AUTOMATIC
    assert_almost_equal(fm.misfit(), 0.617)

    assert mm.type() == "Mww"
    assert_almost_equal(mm.magnitude().value(), 6.366)

    assert_almost_equal(do.depth().value(), 180.5)
    assert_almost_equal(do.latitude().value(), -23.21)
    assert_almost_equal(do.longitude().value(), -66.24)

    nps: NodalPlanes = fm.nodalPlanes()
    np1: NodalPlane = nps.nodalPlane1()
    np2: NodalPlane = nps.nodalPlane2()
    assert_almost_equal(np1.strike().value(), 135.07)
    assert_almost_equal(np1.dip().value(), 21.74)
    assert_almost_equal(np1.rake().value(), 234.4)
    assert_almost_equal(np2.strike().value(), 352.7)
    assert_almost_equal(np2.dip().value(), 72.5)
    assert_almost_equal(np2.rake().value(), -76.9)

    assert fm.momentTensorCount() == 1
    mt: MomentTensor = fm.momentTensor(0)
    assert_almost_equal(mt.clvd(), 0.2206)
    assert_almost_equal(mt.doubleCouple(), 0.7794)
    tensor: Tensor = mt.tensor()
    assert_almost_equal(tensor.Mpp().value(), 2.434e18)
    assert_almost_equal(tensor.Mrp().value(), -3.518e18)
    assert_almost_equal(tensor.Mrr().value(), -2.2605e18)
    assert_almost_equal(tensor.Mrt().value(), 3.085e17)
    assert_almost_equal(tensor.Mtp().value(), -1.381e18)
    assert_almost_equal(tensor.Mtt().value(), -1.735e17)

    assert do.type() == CENTROID
    assert do.quality().usedPhaseCount() == 14
    assert do.quality().usedStationCount() == 14
    assert fm.stationPolarityCount() == 14
    assert mm.stationCount() == 14
    assert_almost_equal(fm.azimuthalGap(), 125.17)

    assert mt.derivedOriginID() == do.publicID()
    assert mm.originID() == do.publicID()


def test_write_scml(golden, tests_dir, tmp_path):
    sample = WPhaseResult.parse_file(tests_dir / "ga2023fsbydd-result.json")
    objects = createObjects(sample, "GA", publicid_slug="test")
    outfile = tmp_path / "new.xml"
    writeSCML(outfile, objects)
    golden.test(tests_dir / "ga2023fsbydd-sc3.xml").check_file(outfile)
