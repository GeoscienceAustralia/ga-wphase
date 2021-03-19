from builtins import str
from builtins import object
from datetime import datetime
import json

SC3_TIME_FORMAT="%Y-%m-%d %H:%M:%S.%f"
class FMItem(object):
    """Represents the focal mechanism resulting from a w-phase run in python-native types.""" 
    def __init__(self):
        self.dip1 = None
        self.dip2 = None
        self.rake1 = None
        self.rake2 = None
        self.str1 = None
        self.str2 = None
        self.lat = None;
        self.lon = None
        self.depth = None
        self.mag = None
        self.mag_type = None
        self.author = ""
        self.tmtp = None
        self.tmtt = None
        self.tmrt = None
        self.tmrr = None
        self.tmrp = None
        self.tmpp = None
        self.usedPhaseCount = None
        self.usedStationCount = None
        self.azimuthalGap = None
        self.originTime = None
        self.centroid = False
        self.overallMisfit = None

def default_log(x):
    pass

class WPhaseParser(object):
    def __init__(self, logger=default_log):
        self.logger = logger
    def read(self, filename=None, json_data=None):
        """Parse a ga-wphase JSON output into a FMItem."""
        log = self.logger
        if filename is not None:
            if json_data is not None:
                raise ValueError('exactly one of filename and json_data may be not None')

            with open(filename) as json_file:
                log("reading file {}".format(filename))
                return self.read(json_data=json.load(json_file))

        else:
            if json_data is None:
                raise ValueError('exactly one of filename and json_data may be not None')

            filename = 'passed json_data'

            if "MomentTensor" not in json_data:
                raise ValueError("MomentTensor is missing: {}".format(filename))

            mt = json_data["MomentTensor"]
            item = FMItem()

            def getter(name, conv, to=None, err=True, dct=None):
                if dct is None:
                    dct = mt

                if name not in dct:
                    msg = '{} is missing in: {}'.format(name, filename)
                    if err: raise ValueError(msg)
                    else: log(msg)

                elif to is None:
                    to = name

                setattr(item, to, conv(dct[name]))

            try:
                getter('dip1',    float)
                getter('dip2',    float)
                getter('rake1',   float)
                getter('rake2',   float)
                getter('str1',    float)
                getter('str2',    float)
                getter('drlat',   float, 'lat')
                getter('drlon',   float, 'lon')
                getter('drdepth', float, 'depth')
                getter('drmag',   float, 'mag')
                getter('drmagt',  str,   'mag_type', False)
                getter('auth',    str,   'author',   False)
                getter('tmtp',    float, err=False)
                getter('tmtt',    float, err=False)
                getter('tmrt',    float, err=False)
                getter('tmrr',    float, err=False)
                getter('tmrp',    float, err=False)
                getter('tmpp',    float, err=False)

                # TODO check if number_of_channels can be really mapped to used phase count
                getter('number_of_channels', int, 'usedPhaseCount',   False, json_data["QualityParams"])
                getter('number_of_stations', int, 'usedStationCount', False, json_data["QualityParams"])
                getter('azimuthal_gap',      float, 'azimuthalGap',     False, json_data["QualityParams"])
                getter(
                    'time',
                    lambda x: datetime.strptime(str(x), SC3_TIME_FORMAT),
                    'originTime',
                    False,
                    json_data["Event"])

                # use OL4 if available, otherwise use OL2, otherwise leave out
                try:
                    item.overallMisfit = float(json_data.get('OL3', json_data.get('OL2'))['misfit']) / 100.
                except Exception:
                    log('could not get misfit')

                item.centroid = "Centroid" in json_data

                log("parsed file %s successfully" % filename)

                return item

            except Exception as e:
                raise ValueError("Failed to load json data: {}".format(e))
