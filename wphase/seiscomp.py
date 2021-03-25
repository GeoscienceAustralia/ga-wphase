"""Methods to convert results to seiscomp formats."""
from contextlib import contextmanager
from seiscomp3 import DataModel as DM, Core, IO
from seiscomp3 import Logging
from wphase.result import FMItem

def datetime_to_seiscomp(dt):
    """Convert a python UTC datetime to a seiscomp Time."""
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
    try:
        wasEnabled = DM.Notifier.IsEnabled()
        DM.Notifier.Enable()
        yield
    finally:
        DM.Notifier.SetEnabled(wasEnabled)


def createObjects(item, agency, evid=None, with_notifiers=False, logging=Logging):
    """Convert an FMItem to seiscomp3.DataModel objects, optionally sending
    Notifier events to messaging.

    :param FMItem item:
        A W-Phase result represented as an FMItem instance.
    :param bool with_notifiers:
        If true, seiscomp3.DataModel.Notifier instances will be created linking
        the FocalMechanism and derived Origin to the triggering event.
    :rtype: dict
    :return:
        A dictionary with keys focalMechanism, derivedOrigin, momentMagnitude
        (mapped to seiscomp3.DataModel objects), and optionally notifiers
        (mapped to a list of Notifiers)."""
    # get current time in UTC
    time = Core.Time.GMT()

    # create creation info
    ci = DM.CreationInfo()
    ci.setCreationTime(time)
    ci.setAgencyID(agency)
    ci.setAuthor(item.author.encode())

    originTime = DM.TimeQuantity(
        datetime_to_seiscomp(item.originTime)
    )

    # fill derived origin
    derivedOrigin = DM.Origin.Create()
    derivedOrigin.setCreationInfo(ci)
    derivedOrigin.setTime(originTime)
    derivedOrigin.setLatitude(DM.RealQuantity(item.lat))
    derivedOrigin.setLongitude(DM.RealQuantity(item.lon))
    derivedOrigin.setDepth(DM.RealQuantity(item.depth))
    derivedOrigin.setEvaluationMode(DM.AUTOMATIC)
    derivedOrigin.setEvaluationStatus(DM.CONFIRMED)

    originQuality = DM.OriginQuality()
    try: originQuality.setUsedPhaseCount(item.usedPhaseCount)
    except Exception: pass

    try: originQuality.setUsedStationCount(item.usedStationCount)
    except Exception: pass

    derivedOrigin.setQuality(originQuality)

    if item.centroid: derivedOrigin.setType(DM.CENTROID)
    else: derivedOrigin.setType(DM.HYPOCENTER)

    # fill magnitude
    try:
        mag = DM.Magnitude.Create()
        mag.setMagnitude(DM.RealQuantity(item.mag))
        mag.setCreationInfo(ci)
        mag.setOriginID(derivedOrigin.publicID())
        mag.setType(item.mag_type.encode())
        mag.setStationCount(item.usedStationCount)
        mag.setMethodID("wphase")
    except Exception as e:
        logging.error('Failed to configure magnitude: {}'.format(e))

    ## Set FocalMechanism
    nodalPlanes = DM.NodalPlanes()
    np1 = DM.NodalPlane()
    np2 = DM.NodalPlane()

    np1.setStrike(DM.RealQuantity(item.str1))
    np1.setDip(DM.RealQuantity(item.dip1))
    np1.setRake(DM.RealQuantity(item.rake1))

    np2.setStrike(DM.RealQuantity(item.str2))
    np2.setDip(DM.RealQuantity(item.dip2))
    np2.setRake(DM.RealQuantity(item.rake2))

    nodalPlanes.setNodalPlane1(np1)
    nodalPlanes.setNodalPlane2(np2)

    fm = DM.FocalMechanism.Create()
    fm.setNodalPlanes(nodalPlanes)
    fm.setCreationInfo(ci)
    fm.setMethodID("wphase")
    fm.setEvaluationMode(DM.AUTOMATIC)
    fm.setMisfit(item.overallMisfit)
    fm.setStationPolarityCount(item.usedPhaseCount)

    try: fm.setAzimuthalGap(item.azimuthalGap)
    except Exception: pass

    # TODO set axis

    # fill tensor
    tensor = DM.Tensor()

    try: tensor.setMtp(DM.RealQuantity(item.tmtp))
    except Exception: pass

    try: tensor.setMtt(DM.RealQuantity(item.tmtt))
    except Exception: pass

    try: tensor.setMrt(DM.RealQuantity(item.tmrt))
    except Exception: pass

    try: tensor.setMrr(DM.RealQuantity(item.tmrr))
    except Exception: pass

    try: tensor.setMrp(DM.RealQuantity(item.tmrp))
    except Exception: pass

    try: tensor.setMpp(DM.RealQuantity(item.tmpp))
    except Exception: pass

    # fill moment tensor object
    mt = DM.MomentTensor.Create()
    mt.setTensor(tensor)
    mt.setCreationInfo(ci)
    mt.setDerivedOriginID(derivedOrigin.publicID())
    mt.setMethodID("wphase")

    try: mt.setMomentMagnitudeID(mag.publicID())
    except Exception: pass

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

def createAndSendObjects(item, connection, logging=Logging, **kwargs):
    """Convert the given FMItem to seiscomp3.DataModel objects, send them over
    the given connection, and return them.

    :param FMItem item: W-Phase result
    :param seiscomp3.Client.Connection connection: seiscomp messaging connection
    :rtype: dict"""
    # create SeiscomP3 objects from focal mechanism item
    with SCNotifier():
        ret = createObjects(item, with_notifiers=True, logging=logging, **kwargs)

    try:
        # serialize objects
        msg = DM.Notifier.GetMessage()

        # forward message to the messaging system
        connection.send(msg)
        logging.info("sent focal mechanism successfully")
    except Exception as e:
        # Log error and continue, since we still created the objs fine
        logging.error('failed send objects to messaging: {}'.format(e))

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
