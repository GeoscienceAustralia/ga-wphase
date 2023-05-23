from functools import partial
from obspy import UTCDateTime
from wphase.psi.model import WPhaseResult
from wphase.seiscomp import createObjects, writeSCML
from seiscomp3.DataModel import AUTOMATIC, CENTROID, FocalMechanism, MomentTensor, NodalPlane, NodalPlanes, Origin, Magnitude, Tensor
import numpy.testing

assert_almost_equal = partial(numpy.testing.assert_almost_equal, decimal=3)

def test_create_objects(tests_dir):
    sample = WPhaseResult.parse_file(tests_dir / "ga2023fsbydd-result.json")
    objects = createObjects(sample, "GA")
    assert set(objects.keys()) == {"focalMechanism", "derivedOrigin", "momentMagnitude"}
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
    assert_almost_equal(fm.misfit(), 0.60476677)

    assert mm.type() == "Mww"
    assert_almost_equal(mm.magnitude().value(), 6.366682797749126)

    assert_almost_equal(do.depth().value(), 170.5)
    assert_almost_equal(do.latitude().value(), -23.213)
    assert_almost_equal(do.longitude().value(), -66.235)

    nps: NodalPlanes = fm.nodalPlanes()
    np1: NodalPlane = nps.nodalPlane1()
    np2: NodalPlane = nps.nodalPlane2()
    assert_almost_equal(np1.strike().value(), 151.0894)
    assert_almost_equal(np1.dip().value(), 22.48)
    assert_almost_equal(np1.rake().value(), 251.977)
    assert_almost_equal(np2.strike().value(), 350.4878)
    assert_almost_equal(np2.dip().value(), 68.6785)
    assert_almost_equal(np2.rake().value(), -82.70401)

    assert fm.momentTensorCount() == 1
    mt: MomentTensor = fm.momentTensor(0)
    assert_almost_equal(mt.clvd(), 0.35437)
    assert_almost_equal(mt.doubleCouple(), 0.64562)
    tensor: Tensor = mt.tensor()
    assert_almost_equal(tensor.Mpp().value(), 3.0810939673170893e+18)
    assert_almost_equal(tensor.Mrp().value(), -3.1716891598093266e+18)
    assert_almost_equal(tensor.Mrr().value(), -2.54973659062008e+18)
    assert_almost_equal(tensor.Mrt().value(), 4.7805636937507776e+17)
    assert_almost_equal(tensor.Mtp().value(), -1.2362758097016248e+18)
    assert_almost_equal(tensor.Mtt().value(), -5.313573766970089e+17)

    assert do.type() == CENTROID
    assert do.quality().usedPhaseCount() == 16
    assert do.quality().usedStationCount() == 16
    assert fm.stationPolarityCount() == 16
    assert mm.stationCount() == 16
    assert_almost_equal(fm.azimuthalGap(), 86.224266)

    assert mt.derivedOriginID() == do.publicID()
    assert mm.originID() == do.publicID()
