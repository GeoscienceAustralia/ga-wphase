"""Data models for W-Phase calculation input and outputs.

There is no good naming convention here - unfortunately, the JSON output must
follow a fixed schema to work with other GA systems, so we're kinda stuck with
it."""


from typing import Callable, List, Literal, Mapping, Optional, OrderedDict, Tuple, TYPE_CHECKING, Union

import numpy as np
from obspy import UTCDateTime as _UTCDateTime
from pydantic import BaseModel, Field

# Sadly we have to jump through some hoops to cleanly add pydantic validators
# to third-party types if we want it to type-check correctly:
if TYPE_CHECKING:
    UTCDateTime = _UTCDateTime
    Array = np.ndarray
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
        json_encoders: Mapping[type, Callable] = {
            _UTCDateTime: lambda t: str(t).replace("T", " ").replace("Z", ""),
            np.ndarray: np.ndarray.tolist,
            complex: lambda z: (z.imag, z.real)
        }
        # since we have some non-serialized ndarray fields:
        arbitrary_types_allowed = True
        extra = "forbid"


class PreliminaryMagnitudeFit(Data):
    magnitude: float
    unclamped_magnitude: float
    scalar_moment: float
    timescale: float
    regularization: float
    strike: float
    average_amplitude: float
    anisotropy: float
    eccentricity: float
    corrected_amplitudes: np.ndarray
    azimuths: np.ndarray
    trids: List[str]


class OL1Result(Data):
    """Result of the preliminary magnitude fit, which is used as a first guess
    for the actual W-Phase inversion."""

    magnitude: float
    """Mw magnitude estimated using station azimuths, distances and p2p amplitude."""
    nstations: int
    """Number of stations used in fit"""
    used_traces: List[str]
    """List of waveform stream IDs of channels used in fit"""

    preliminary_calc_details: PreliminaryMagnitudeFit = Field(exclude=True)


class OL2Result(Data):
    """A result of the first inversion for the deviatoric moment tensor.

    This model includes serialized fields (which are used to produce the wphase
    output JSON) and non-serialized fields (which are for internal use only and
    contain heavy-duty data like the synthetic waveforms)."""

    Mrr: float
    """R-R component of moment tensor in Newton-metres"""
    Mtt: float
    """Θ-Θ component of moment tensor in Newton-metres"""
    Mpp: float
    """Φ-Φ component of moment tensor in Newton-metres"""
    Mrt: float
    """R-Θ component of moment tensor in Newton-metres"""
    Mrp: float
    """R-Φ component of moment tensor in Newton-metres"""
    Mtp: float
    """Θ-Φ component of moment tensor in Newton-metres"""
    misfit: float
    """Solution misfit (percentage L2 difference between observed and synthetic
    displacements)"""
    m0: float
    """Scalar moment in Newton-metres"""
    magnitude: float
    """Mww magnitude"""
    depth: float
    """Depth in km used in inversion"""
    time_delay: float
    """Time delay in seconds"""
    used_traces: List[str]
    """List of waveform stream IDs of channels used in inversion"""
    trace_lengths: OrderedDict[str, int]
    """Length of each waveform trace, keyed on waveform stream ID"""
    periods: Tuple[float, float]
    """Periods of the bandpass filter used (low, high; in seconds)"""

    moment_tensor: np.ndarray = Field(exclude=True)
    """Moment tensor in Newton-metres as a 6-element array"""
    observed_displacements: np.ndarray = Field(exclude=True)
    """Observed displacements (after filtering)"""
    synthetic_displacements: np.ndarray = Field(exclude=True)
    """Synthetic displacements for this solution"""


class OL3Result(OL2Result):
    """Like OL2, but solved for not only the moment tensor but also the
    centroid."""

    centroid: Tuple[float, float, float]
    """Optimal centroid location as tuple (latitude (deg), longitude (deg), depth (km))"""

    grid_search_candidates: Optional[List[Tuple[float, float, float]]] = Field(exclude=True)
    """List of inputs to core_inversion that were used in the grid search.
    Elements are (lat, lon, depths) tuples."""

    grid_search_results: Optional[List[Tuple[np.ndarray, float]]] = Field(exclude=True)
    """List of outputs from the inversion for each point in the grid search.
    elements are (MT, misfit) tuples."""


class Quality(Data):
    """Basic quality parameters of a W-Phase solution."""

    azimuthal_gap: float
    number_of_stations: float
    number_of_channels: float


class Event(Data):
    """The seismic event parameters used as the first guess for an inversion."""

    id: Optional[str] = None
    depth: float
    latitude: float
    longitude: float
    time: UTCDateTime


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


class TimeDelayMisfits(Data):
    array: List[float]
    min: int


class Timing(Data):
    ol1: Optional[float] = None
    ol2: Optional[float] = None
    ol3: Optional[float] = None
    new: Optional[float] = None


class WPhaseResult(Data):
    """This defines the schema of ga-wphase's final JSON output."""

    Event: Event

    OL1: Optional[OL1Result] = None
    OL2: Optional[OL2Result] = None
    OL3: Optional[OL3Result] = None
    new: Optional[OL3Result] = None
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
    timing: Timing = Timing()

    DataSource: Optional[str] = None
    HostName: Optional[str] = None

    def add_warning(self, warning):
        self.Warnings.append(str(warning))


class ChannelMetadata(Data):
    """Metadata for a single waveform channel."""
    azimuth: float
    dip: float
    elevation: float
    gain: float
    latitude: float
    """Latitude of station in degrees"""
    longitude: float
    """Longitude of station in degrees"""
    sensitivity: float
    """Seismometer sensitivity"""
    transfer_function: Union[Literal["A"], Literal["B"], None]
    """Type (A or B) of the transfer function.

    Type A uses frequency units of rads/second, while Type B uses Hz."""
    zeros: List[complex]
    """Zeroes of the transfer function"""
    poles: List[complex]
    """Poles of the transfer function"""
    sampling_rate: float
    """Samples per second"""
