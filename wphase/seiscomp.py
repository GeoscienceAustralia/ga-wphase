"""Methods to convert results to seiscomp formats."""
from contextlib import contextmanager
import logging


from obspy.core import UTCDateTime
from seiscomp3 import DataModel as DM, Core, IO
from seiscomp3 import Logging
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
            Logging.debug(b"Detected SWIG char* type as bytes")
            _charstar_is_bytes = True
        except TypeError:
            Logging.debug(u"Detecting SWIG char* type as unicode")
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
    return Core.Time(dt.year,
                     dt.month,
                     dt.day,
                     dt.hour,
                     dt.minute,
                     dt.second,
                     dt.microsecond)


@contextmanager
def SCNotifier():
    """seiscomp entities created (and associated with each other) inside this
    context manager will be logged to the Notifier."""
    wasEnabled = DM.Notifier.IsEnabled()
    try:
        DM.Notifier.Enable()
        yield
    finally:
        DM.Notifier.SetEnabled(wasEnabled)


def createObjects(item: WPhaseResult, agency, evid=None, with_notifiers=False):
    """Convert a WPhaseResult to seiscomp3.DataModel objects, optionally sending
    Notifier events to messaging.

    :param WPhaseResult item:
        Result of a W-Phase inversion.
    :param bool with_notifiers:
        If true, seiscomp3.DataModel.Notifier instances will be created linking
        the FocalMechanism and derived Origin to the triggering event.
    :rtype: dict
    :return:
        A dictionary with keys focalMechanism, derivedOrigin, momentMagnitude
        (mapped to seiscomp3.DataModel objects), and optionally notifiers
        (mapped to a list of Notifiers)."""
    mt = item.MomentTensor
    preferredOL = item.OL3 or item.OL2
    if not (mt and preferredOL):
        raise ValueError("W-Phase result contains no MomentTensor to convert/send!")

    # get current time in UTC
    time = Core.Time.GMT()

    # create creation info
    ci = DM.CreationInfo()
    ci.setCreationTime(time)
    ci.setAgencyID(agency)
    ci.setAuthor(charstar(mt.auth))

    originTime = DM.TimeQuantity(datetime_to_seiscomp(item.Event.time))

    # fill derived origin
    derivedOrigin = DM.Origin.Create()
    derivedOrigin.setCreationInfo(ci)
    derivedOrigin.setTime(originTime)
    derivedOrigin.setLatitude(DM.RealQuantity(mt.drlat))
    derivedOrigin.setLongitude(DM.RealQuantity(mt.drlon))
    derivedOrigin.setDepth(DM.RealQuantity(mt.drdepth))
    derivedOrigin.setEvaluationMode(DM.AUTOMATIC)
    derivedOrigin.setEvaluationStatus(DM.CONFIRMED)

    originQuality = DM.OriginQuality()
    if item.QualityParams is not None:
        try:
            originQuality.setUsedPhaseCount(item.QualityParams.number_of_channels)
        except Exception:
            pass

        try:
            originQuality.setUsedStationCount(item.QualityParams.number_of_stations)
        except Exception:
            pass

    derivedOrigin.setQuality(originQuality)

    if item.Centroid:
        derivedOrigin.setType(DM.CENTROID)
    else:
        derivedOrigin.setType(DM.HYPOCENTER)

    # fill magnitude
    try:
        mag = DM.Magnitude.Create()
        mag.setMagnitude(DM.RealQuantity(mt.drmag))
        mag.setCreationInfo(ci)
        mag.setOriginID(derivedOrigin.publicID())
        mag.setType(charstar(mt.drmagt))
        if item.QualityParams:
            mag.setStationCount(item.QualityParams.number_of_stations)
        mag.setMethodID("wphase")
    except Exception as e:
        logger.error("Failed to configure magnitude: {}".format(e))

    ## Set FocalMechanism
    nodalPlanes = DM.NodalPlanes()
    np1 = DM.NodalPlane()
    np2 = DM.NodalPlane()

    np1.setStrike(DM.RealQuantity(mt.str1))
    np1.setDip(DM.RealQuantity(mt.dip1))
    np1.setRake(DM.RealQuantity(mt.rake1))

    np2.setStrike(DM.RealQuantity(mt.str2))
    np2.setDip(DM.RealQuantity(mt.dip2))
    np2.setRake(DM.RealQuantity(mt.rake2))

    nodalPlanes.setNodalPlane1(np1)
    nodalPlanes.setNodalPlane2(np2)

    fm = DM.FocalMechanism.Create()
    fm.setNodalPlanes(nodalPlanes)
    fm.setCreationInfo(ci)
    fm.setMethodID("wphase")
    fm.setEvaluationMode(DM.AUTOMATIC)

    misfit = preferredOL.misfit / 100.0
    fm.setMisfit(misfit)
    if item.QualityParams:
        fm.setStationPolarityCount(item.QualityParams.number_of_channels)
        try:
            fm.setAzimuthalGap(item.QualityParams.azimuthal_gap)
        except Exception:
            pass

    # TODO set axis

    # fill tensor
    tensor = DM.Tensor()

    try:
        tensor.setMtp(DM.RealQuantity(fm.tmtp))
    except Exception:
        pass

    try:
        tensor.setMtt(DM.RealQuantity(fm.tmtt))
    except Exception:
        pass

    try:
        tensor.setMrt(DM.RealQuantity(fm.tmrt))
    except Exception:
        pass

    try:
        tensor.setMrr(DM.RealQuantity(fm.tmrr))
    except Exception:
        pass

    try:
        tensor.setMrp(DM.RealQuantity(fm.tmrp))
    except Exception:
        pass

    try:
        tensor.setMpp(DM.RealQuantity(fm.tmpp))
    except Exception:
        pass

    # fill moment tensor object
    mt = DM.MomentTensor.Create()
    mt.setTensor(tensor)
    mt.setCreationInfo(ci)
    mt.setDerivedOriginID(derivedOrigin.publicID())
    mt.setMethodID("wphase")

    try:
        mt.setClvd(fm.clvd)
    except Exception:
        pass

    try:
        mt.setDoubleCouple(fm.dc)
    except Exception:
        pass

    try:
        mt.setMomentMagnitudeID(mag.publicID())
    except Exception:
        pass

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

    # Adding these seems to *immediately* queue up the notifications; so we
    # must do this *after* notifying the origin/FM so that scmaster
    # processes everything in the correct order.
    derivedOrigin.add(mag)
    fm.add(mt)

    ret = {
        "derivedOrigin": derivedOrigin,
        "momentMagnitude": mag,
        "focalMechanism": fm,
    }
    if with_notifiers:
        ret["notifiers"] = notifiers
    return ret


def createAndSendObjects(item, connection, **kwargs):
    """Convert the given FMItem to seiscomp3.DataModel objects, send them over
    the given connection, and return them.

    :param FMItem item: W-Phase result
    :param seiscomp3.Client.Connection connection: seiscomp messaging connection
    :rtype: dict"""
    # create SeiscomP3 objects from focal mechanism item
    with SCNotifier():
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


def writeSCML(filename, objects):
    """Given seiscomp3.DataModel objects (as produced by createObjects),
    write them to an XML file.

    :param str filename: path to output file
    :param objects dict: same as return type of createObjects"""
    # create SeisComP3 XML Archive used to serialize objects
    ar = IO.XMLArchive()

    # enable formatted output
    ar.setFormattedOutput(True)

    # try to create the output file
    ar.create(filename)

    # Serialize the objects
    for x in objects.values():
        if isinstance(x, DM.PublicObject):
            ar.writeObject(x)

    ar.close()
    return True
