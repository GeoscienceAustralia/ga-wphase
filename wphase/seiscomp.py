"""Methods to convert results to seiscomp formats."""
from __future__ import annotations
from contextlib import contextmanager
import logging
from pathlib import Path
from typing import List, Optional, Union

from obspy.core import UTCDateTime
from seiscomp import datamodel as DM, core, io
import seiscomp.logging
from wphase.psi.model import WPhaseResult

logger = logging.getLogger(__name__)

_charstar_is_bytes = None

def charstar(string):
    """Convert a string (unicode in python3, bytes in python2) to a char*
    usable as an argument to the seiscomp SWIG API.

    Depending on what version of seiscomp and python we're using, and whether
    seiscomp's SWIG bindings were generated with SWIG_PYTHON_2_UNICODE or not,
    the correct type to feed to C++ string/char* arguments can vary.
    Unfortunately, the seiscomp3 backwards-compat python wrapper doesn't
    compensate for this.

    I couldn't find a simple way to introspect the correct python type, so the
    first time this method is called it attempts to log a message to find out."""
    global _charstar_is_bytes
    if _charstar_is_bytes is None:
        # first time we've been called - we need to detect.
        try:
            seiscomp.logging.debug(b"Detected SWIG char* type as bytes")
            _charstar_is_bytes = True
        except TypeError:
            seiscomp.logging.debug(u"Detecting SWIG char* type as unicode")
            _charstar_is_bytes = False
    if _charstar_is_bytes:
        if isinstance(string, bytes):
            return string
        else:
            return string.encode('utf-8')
    else:
        if isinstance(string, bytes):
            return string.decode('utf-8')
        else:
            return string

def datetime_to_seiscomp(dt):
    """Convert a python or obspy UTC datetime to a seiscomp Time."""
    if isinstance(dt, UTCDateTime):
        dt = dt.datetime
    return core.Time(dt.year,
                     dt.month,
                     dt.day,
                     dt.hour,
                     dt.minute,
                     dt.second,
                     dt.microsecond)


@contextmanager
def SCNotifierEnabled():
    """seiscomp entities created (and associated with each other) inside this
    context manager will be logged to the Notifier."""
    wasEnabled = DM.Notifier.IsEnabled()
    try:
        DM.Notifier.Enable()
        yield
    finally:
        DM.Notifier.SetEnabled(wasEnabled)


try:
    from typing import TypedDict
except ImportError:
    pass
else:

    class SeiscompResults(TypedDict):
        momentMagnitude: DM.Magnitude
        derivedOrigin: DM.Origin
        focalMechanism: DM.FocalMechanism
        notifiers: Optional[List[DM.Notifier]]


def createObjects(
    item: WPhaseResult,
    agency,
    evid=None,
    with_notifiers=False,
    publicid_slug: Optional[str] = None,
    triggering_origin_id: Optional[str] = None,
) -> SeiscompResults:
    """Convert a WPhaseResult to seiscomp.datamodel objects, optionally sending
    Notifier events to messaging.

    :param WPhaseResult item:
        Result of a W-Phase inversion.
    :param bool with_notifiers:
        If true, seiscomp.datamodel.Notifier instances will be created linking
        the FocalMechanism and derived Origin to the triggering event.
    :param str publicid_slug:
        If provided, publicIDs are built from templates using this string
        rather than dynamically based on time. (Used for deterministic outputs
        in tests.)
    :rtype: dict
    :return:
        A dictionary with keys focalMechanism, derivedOrigin, momentMagnitude
        (mapped to seiscomp.datamodel objects), and optionally notifiers
        (mapped to a list of Notifiers)."""
    mtresult = item.MomentTensor
    preferredOL = item.OL3 or item.OL2
    if not (mtresult and preferredOL):
        raise ValueError("W-Phase result contains no MomentTensor to convert/send!")

    if item.CreationTime is None:
        # default to current time in UTC
        time = core.Time.GMT()
    else:
        time = datetime_to_seiscomp(item.CreationTime)

    # create creation info
    ci = DM.CreationInfo()
    ci.setCreationTime(time)
    ci.setAgencyID(agency)
    ci.setAuthor(charstar(mtresult.auth))

    originTime = DM.TimeQuantity(datetime_to_seiscomp(item.Event.time))

    # fill derived origin
    derivedOrigin = DM.Origin.Create()
    if publicid_slug is not None:
        derivedOrigin.setPublicID(f"Origin/{publicid_slug}")
    derivedOrigin.setCreationInfo(ci)
    derivedOrigin.setTime(originTime)
    derivedOrigin.setLatitude(DM.RealQuantity(mtresult.drlat))
    derivedOrigin.setLongitude(DM.RealQuantity(mtresult.drlon))
    derivedOrigin.setDepth(DM.RealQuantity(mtresult.drdepth))
    derivedOrigin.setEvaluationMode(DM.AUTOMATIC)
    derivedOrigin.setEvaluationStatus(DM.CONFIRMED)

    originQuality = DM.OriginQuality()
    if item.QualityParams is not None:
        try:
            originQuality.setUsedPhaseCount(int(item.QualityParams.number_of_channels))
        except Exception:
            pass

        try:
            originQuality.setUsedStationCount(
                int(item.QualityParams.number_of_stations)
            )
        except Exception:
            pass

    derivedOrigin.setQuality(originQuality)
    derivedOrigin.setMethodID("wphase")

    if item.Centroid:
        derivedOrigin.setType(DM.CENTROID)
    else:
        derivedOrigin.setType(DM.HYPOCENTER)

    # fill magnitude
    mag = DM.Magnitude.Create()
    if publicid_slug is not None:
        mag.setPublicID(f"Magnitude/{publicid_slug}")
    try:
        mag.setMagnitude(DM.RealQuantity(mtresult.drmag))
        mag.setCreationInfo(ci)
        mag.setOriginID(derivedOrigin.publicID())
        mag.setType(charstar(mtresult.drmagt))
        if item.QualityParams:
            mag.setStationCount(int(item.QualityParams.number_of_stations))
        mag.setMethodID("wphase")
    except Exception as e:
        logger.error("Failed to configure magnitude: {}".format(e))

    ## Set FocalMechanism
    nodalPlanes = DM.NodalPlanes()
    np1 = DM.NodalPlane()
    np2 = DM.NodalPlane()

    np1.setStrike(DM.RealQuantity(mtresult.str1))
    np1.setDip(DM.RealQuantity(mtresult.dip1))
    np1.setRake(DM.RealQuantity(mtresult.rake1))

    np2.setStrike(DM.RealQuantity(mtresult.str2))
    np2.setDip(DM.RealQuantity(mtresult.dip2))
    np2.setRake(DM.RealQuantity(mtresult.rake2))

    nodalPlanes.setNodalPlane1(np1)
    nodalPlanes.setNodalPlane2(np2)

    fm: DM.FocalMechanism = DM.FocalMechanism.Create()
    if publicid_slug is not None:
        fm.setPublicID(f"FocalMechanism/{publicid_slug}")
    if triggering_origin_id is not None:
        fm.setTriggeringOriginID(triggering_origin_id)
    fm.setNodalPlanes(nodalPlanes)
    fm.setCreationInfo(ci)
    fm.setMethodID("wphase")
    fm.setEvaluationMode(DM.AUTOMATIC)

    misfit = preferredOL.misfit / 100.0
    fm.setMisfit(misfit)
    if item.QualityParams:
        fm.setStationPolarityCount(int(item.QualityParams.number_of_channels))
        try:
            fm.setAzimuthalGap(item.QualityParams.azimuthal_gap)
        except Exception:
            pass

    # TODO set axis

    # fill tensor
    tensor = DM.Tensor()

    try:
        tensor.setMtp(DM.RealQuantity(mtresult.tmtp))
    except Exception:
        pass

    try:
        tensor.setMtt(DM.RealQuantity(mtresult.tmtt))
    except Exception:
        pass

    try:
        tensor.setMrt(DM.RealQuantity(mtresult.tmrt))
    except Exception:
        pass

    try:
        tensor.setMrr(DM.RealQuantity(mtresult.tmrr))
    except Exception:
        pass

    try:
        tensor.setMrp(DM.RealQuantity(mtresult.tmrp))
    except Exception:
        pass

    try:
        tensor.setMpp(DM.RealQuantity(mtresult.tmpp))
    except Exception:
        pass

    # fill moment tensor object
    mt = DM.MomentTensor.Create()
    if publicid_slug is not None:
        mt.setPublicID(f"MomentTensor/{publicid_slug}")
    mt.setTensor(tensor)
    mt.setCreationInfo(ci)
    mt.setDerivedOriginID(derivedOrigin.publicID())
    mt.setMethodID("wphase")

    try:
        mt.setClvd(mtresult.clvd)
    except Exception:
        pass

    try:
        mt.setDoubleCouple(mtresult.dc)
    except Exception:
        pass

    try:
        mt.setMomentMagnitudeID(mag.publicID())
    except Exception:
        pass

    for tr in preferredOL.used_traces:
        try:
            n, s, l, c = tr.split(".")
            wfid = DM.WaveformStreamID(n, s, l, c, "")
            comp = DM.MomentTensorComponentContribution()
            comp.setComponent(0)
            comp.setActive(True)
            comp.setPhaseCode("W")
            if (misfit_pc := preferredOL.trace_misfits.get(tr)) is not None:
                comp.setMisfit(misfit_pc/100)
            comp.setWeight(1)
            cont = DM.MomentTensorStationContribution.Create()
            if publicid_slug is not None:
                cont.setPublicID(f"MomentTensorStationContribution/{publicid_slug}.{n}.{s}.{l}.{c}")
            cont.setActive(True)
            cont.setWeight(1)
            cont.setWaveformID(wfid)
            cont.add(comp)
        except Exception:
            logger.exception(
                "Error constructing momentTensorStationContribution for "
                "%s, omitting this station",
                tr
            )
            assert False
        else:
            mt.add(cont)


    # Since we don't want to overwrite the event data itself, but we do
    # want to explicitly associate our data to the correct event, we have
    # to manually create Notifiers for these associations.
    oRef = DM.OriginReference()
    oRef.setOriginID(derivedOrigin.publicID())

    fmRef = DM.FocalMechanismReference()
    fmRef.setFocalMechanismID(fm.publicID())
    if with_notifiers and evid:
        notifiers = [
            # TODO: are these two actually valid? Origins/FMs are associated to events
            # by references, but are not children thereof!
            DM.Notifier.Create(evid, DM.OP_ADD, derivedOrigin),
            DM.Notifier.Create(evid, DM.OP_ADD, fm),
            # I *think* we only actually need these two.
            DM.Notifier.Create(evid, DM.OP_ADD, oRef),
            DM.Notifier.Create(evid, DM.OP_ADD, fmRef),
        ]
    else:
        notifiers = None

    # Adding these seems to *immediately* queue up the notifications; so we
    # must do this *after* notifying the origin/FM so that scmaster
    # processes everything in the correct order.
    derivedOrigin.add(mag)
    fm.add(mt)

    return {
        "derivedOrigin": derivedOrigin,
        "momentMagnitude": mag,
        "focalMechanism": fm,
        "notifiers": notifiers,
    }


def createAndSendObjects(item, connection, **kwargs) -> SeiscompResults:
    """Convert the given FMItem to seiscomp.datamodel objects, send them over
    the given connection, and return them.

    :param FMItem item: W-Phase result
    :param seiscomp.client.Connection connection: seiscomp messaging connection
    :rtype: dict"""
    # create seiscomp objects from focal mechanism item
    with SCNotifierEnabled():
        ret = createObjects(item, with_notifiers=True, **kwargs)

    try:
        # serialize objects
        msg = DM.Notifier.GetMessage()

        # forward message to the messaging system
        connection.send(msg)
        logger.info("sent focal mechanism successfully")
    except Exception as e:
        # Log error and continue, since we still created the objs fine
        logger.error("failed send objects to messaging: {}".format(e))

    return ret


def writeSCML(filename: Union[str, Path], objects: SeiscompResults):
    """Given seiscomp.datamodel objects (as produced by createObjects),
    write them to an XML file.

    :param str filename: path to output file
    :param objects dict: same as return type of createObjects"""
    # create seiscomp XML Archive used to serialize objects
    ar = io.XMLArchive()

    # enable formatted output
    ar.setFormattedOutput(True)

    # try to create the output file
    ar.create(str(filename))

    ep = DM.EventParameters()
    ep.add(objects["derivedOrigin"])
    ep.add(objects["focalMechanism"])
    ar.writeObject(ep)
    ar.close()
