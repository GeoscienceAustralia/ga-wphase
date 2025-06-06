"""Data models for W-Phase calculation input and outputs.

There is no good naming convention here - unfortunately, the JSON output must
follow a fixed schema to work with other GA systems, so we're kinda stuck with
it."""


from typing import Dict, List, Optional, OrderedDict, Tuple, TYPE_CHECKING

import numpy as np
from obspy import UTCDateTime as _UTCDateTime
from pydantic import BaseModel, Field

# Sadly we have to jump through some hoops to cleanly add pydantic validators
# to third-party types if we want it to type-check correctly:
if TYPE_CHECKING:
    from obspy import UTCDateTime
    from numpy import ndarray as Array
else:
    class UTCDateTime(_UTCDateTime):
        @classmethod
        def __get_validators__(cls):
            yield lambda v: _UTCDateTime(v)

    class Array(np.ndarray):
        @classmethod
        def __get_validators__(cls):
            yield np.asarray


class Data(BaseModel):
    class Config:
        json_encoders = {
            # We keep this datetime format for backwards compatibility:
            _UTCDateTime: lambda t: str(t).replace("T", " ").replace("Z", ""),
            np.ndarray: np.ndarray.tolist,
        }
        # since we have some non-serialized ndarray fields:
        arbitrary_types_allowed = True
        extra = "forbid"


class OL1Result(Data):
    """Result of the preliminary magnitude fit, which is used as a first guess
    for the actual W-Phase inversion."""

    magnitude: float
    nstations: int
    used_traces: List[str]

    preliminary_calc_details: Optional[dict] = Field(exclude=True)


class ScreeningStage(Data):
    name: str
    station_results: Dict[str, bool]
    info: str = ""
    mt: Optional[List[float]] = None


class OL2Result(Data):
    """A result of the first inversion for the deviatoric moment tensor.

    This model includes serialized fields (which are used to produce the wphase
    output JSON) and non-serialized fields (which are for internal use only and
    contain heavy-duty data like the synthetic waveforms)."""

    Mrr: float
    Mtt: float
    Mpp: float
    Mrt: float
    Mrp: float
    Mtp: float
    misfit: float
    m0: float
    magnitude: float
    depth: float
    time_delay: float
    used_traces: List[str]
    trace_lengths: OrderedDict[str, int]
    trace_misfits: OrderedDict[str, float]
    trace_time_windows: Optional[OrderedDict[str, Tuple[float, float]]]

    moment_tensor: Optional[np.ndarray] = Field(exclude=True)
    observed_displacements: Optional[np.ndarray] = Field(exclude=True)
    synthetic_displacements: Optional[np.ndarray] = Field(exclude=True)


class OL3Result(OL2Result):
    """Like OL2, but solved for not only the moment tensor but also the
    centroid."""

    centroid: Tuple[float, float, float]
    centroid_time: UTCDateTime

    grid_search_candidates: Optional[List[Tuple[float, float, float]]] = Field(exclude=True)
    """List of inputs to core_inversion that were used in the grid search.
    Elements are (lat, lon, depths) tuples."""

    grid_search_results: Optional[List[Tuple[np.ndarray, float]]] = Field(exclude=True)
    """List of outputs from the inversion for each point in the grid search.
    elements are (MT, misfit) tuples."""


class Quality(Data):
    """Basic quality parameters of a W-Phase solution."""

    azimuthal_gap: float
    number_of_stations: int
    number_of_channels: int


class Event(Data):
    """The seismic event parameters used as the first guess for an inversion."""

    id: Optional[str] = None
    depth: float
    latitude: float
    longitude: float
    time: UTCDateTime
    creation_time: Optional[UTCDateTime] = None


class AntelopeMomentTensor(Data):
    """A seismic moment tensor expressed in the format used by Antelope,
    apparently."""

    # The upper-triangular components of the moment tensor
    tmpp: float
    tmrp: float
    tmrr: float
    tmrt: float
    tmtp: float
    tmtt: float

    scm: float
    """Seismic moment"""
    drmag: float
    """Magnitude value"""
    drmagt: str
    """Magnitude type"""

    drlat: float
    drlon: float
    drdepth: float

    str1: float
    dip1: float
    rake1: float
    str2: float
    dip2: float
    rake2: float

    dc: Optional[float] = None
    """Double-couple component (from 0 to 1)"""
    clvd: Optional[float] = None
    """Compensated linear vector dipole component (from 0 to 1)"""

    auth: str


class CentroidLocation(Data):
    """The earthquake's location as estimated by the centroid of the waveform
    inversion."""

    depth: float
    latitude: float
    longitude: float
    time: UTCDateTime


class TimeDelayMisfits(Data):
    array: List[float]
    min: int


class WPhaseResult(Data):
    """This defines the schema of ga-wphase's final JSON output."""

    Event: Event

    OL1: Optional[OL1Result] = None
    OL2: Optional[OL2Result] = None
    OL3: Optional[OL3Result] = None
    QualityParams: Optional[Quality] = None
    MomentTensor: Optional[AntelopeMomentTensor] = None  # "Preferred solution"
    Centroid: Optional[CentroidLocation] = None
    misfits: Optional[TimeDelayMisfits] = None

    WphaseResultPlots: Optional[List[Tuple[Tuple[int, int], str]]] = None
    """List of waveform plots produced. Each element is a pair, with the first
    element giving the range of trace indices included in the plot and the
    second element giving the filename."""

    Warnings: List[str] = []
    Error: Optional[str] = None
    StackTrace: Optional[str] = None
    WPInvProfile: Optional[str] = None

    DataSource: Optional[str] = None
    HostName: Optional[str] = None
    CreationTime: Optional[UTCDateTime] = None

    ScreeningStages: List[ScreeningStage] = []

    def add_warning(self, warning):
        self.Warnings.append(str(warning))
