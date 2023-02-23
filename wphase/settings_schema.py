# -*- coding: utf-8 -*-

"""Configuration settings for :mod:`wphase`. These can be adjusted by setting the
corresponding environment variables (prefixed with WPHASE_)."""

import os
from typing import List, Literal, Optional, Tuple, Union

from pydantic import BaseSettings, Field, validator

# We seem to get *better* performance by disabling multithreading in openblas,
# even when we don't use any other parallelism.  Not really sure on the cause,
# but it's true; so here we are:
os.environ['OPENBLAS_NUM_THREADS'] = '1'

class WPhaseSettings(BaseSettings):
    class Config:
        env_prefix = 'WPHASE_'

    ### Paths
    OUTPUT_DIR: str = '/tmp/wphase-output'
    """Default output directory"""

    # This *must* be supplied as environment variable. We define it using Field
    # to avoid spurious type-checking errors:
    GREENS_FUNCTIONS: str = Field(str)
    """Path to greens functions, either a directory or a hdf5 archive"""

    ### Metadata
    AUTHORITY: str = "GA W-phase"
    """Authority name to embed in output"""

    HOST_NAME: str = "localhost"
    """Hostname to embed in output"""

    @validator("AUTHORITY", "HOST_NAME")
    def disallow_empty_string(cls, value):
        if value == "":
            raise ValueError("String must not be empty")
        return value

    ### Model parameters and quality thresholds
    WAVEFORM_COMPONENTS: Optional[List[str]] = ["Z"]
    """List of waveform components to use in inversion. If empty, all components will be used."""

    WORKER_COUNT: Optional[int] = None
    """Number of worker processes to use in parallel computations"""

    BANDPASS_IMPLEMENTATION: Union[Literal["scipy"], Literal["fortran"]] = "scipy"
    """Name of the bandpass filter implementation to use (scipy or fortran)"""

    PROFILE: bool = False
    """Whether performance profiling should be enabled"""

    AMPLITUDE_REGULARIZATION: float = 0.2
    """Strength of regularization in preliminary magnitude calculation.

    Higher values result in less weight being given to azimuthal variation - in
    the limit at infinity, the azimuths are completely ignored.

    0.2 seems to work well enough to prevent issues without significantly
    suppressing anisotropy in datasets with good azimuthal coverage."""

    MINIMUM_PRELIMINARY_MAGNITUDE: float = 6.5
    """Floor to clamp preliminary magnitude to"""

    RESPONSE_MISFIT_TOL: float = 5. # Percent
    """Misfit tolerance in percent.

    A channel is rejected if the fitting of the instrument response is greater
    than this.  See Kanamori2008 for details on the fitting of the instrument
    response."""

    WPHASE_CUTOFF: float = 15
    """Time cutoff factor: if we call this k, then the selected time window for
    each trace will end at t_p + k*d where t_p is the theoretical p-wave
    arrival time and d is the hypocentral distance in degrees."""

    MEDIAN_REJECTION_COEFF: Tuple[float, float] = (0.1, 3)
    """Traces with peak-to-peak amplitude outside of the range
    [median_rejection_coeff[0] * m, median_rejection_coeff[1] * m], where m
    is the median peak-to-peak amplitude over all stations, will be rejected"""

    MINIMUM_STATIONS: int = 4
    """Minimum number of stations required to begin the inversion."""

    MAXIMUM_AZIMUTHAL_GAP: float = 360
    """Maximum azimuthal gap allowed to begin the inversion (in degrees)

    By default this is set to 360 and thus does nothing."""

    MAXIMUM_TIME_DELAY: float = 200.
    """Maximum time delay t_d in seconds. We'll search for the optimal t_d up
    to this value."""

    MISFIT_TOL_SEQUENCE: List[float] = [300, 200, 100]
    """List of thresholds used in outlier removal. Channels where the misfits
    (100*sqrt(sum(synthetic-observed)^2 / sum(synthetic)^2)) are greater than
    the first element of this sequence will be rejected. If any stations
    remain, then the inversion will be repeated with these stations, and those
    with misfits greater than the second element will be rejected, and so
    on."""

    OPTIMIZATION_TOLERANCE: float = 0.5
    """`tol` parameter to pass to scipy.optimize.minimize when performing
    outlier removal. Values above 1 produce bad results. Smaller values produce
    better results but increase runtime."""

    MINIMUM_FITTING_CHANNELS: float = 10
    """Minimum number of well-fitting channels required for inversion to pass
    quality check"""



    ### Output products:

    N_TRACES_PER_RESULT_PLOT: int = 6
    """Number of traces to show in each waveform plot image."""

    OUTPUT_FILE_NAME = 'wphase_output.json'
    """The name of the JSON file to contain the output from wphase calculations."""

    BEACHBALL_PREFIX = "beachball"
    """Prefix for the names of the images containing the beachball."""

    STATION_DISTRIBUTION_PREFIX = "station_distribution"
    """Prefix for the station distribution plot."""

    GRID_SEARCH_PREFIX = "grid_search"
    """Prefix for the grid search plot."""

    STATION_TRACES_PREFIX = "raw_traces"
    """Prefix for the names of the images containing the traces."""

    RESULTS_TRACES_PREFIX = "wphase_traces"
    """Prefix for the names of the images containing the wphase traces."""

    MISFIT_PREFIX = "wphase_misfit"
    """Prefix for the names of the images containing the wphase misfit info."""

    PRELIM_FIT_PREFIX = "preliminary_magnitude_fit"
    """Prefix for the names of the images containing the preliminary magnitude fit"""
