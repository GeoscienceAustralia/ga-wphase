# -*- coding: utf-8 -*-

"""
.. module:: settings

Configuration settings for :mod:`wphase`. At present, the settings are declared
directly as python variables, but those that should be changed by users can
also be overridden with environment variables."""

from __future__ import absolute_import

import os
import sys
import multiprocessing as mp

# Set the number of threads used by OpenBLAS.
os.environ['OPENBLAS_NUM_THREADS'] = '1'

#: Directory containing the green's functions.
GREEN_DIR = os.environ['WPHASE_GREENS_FUNCTIONS']

DEFAULT_OUTPUT_DIR = os.environ.get('WPHASE_OUTPUT_DIR', '/tmp/wphase-output')

#: The number of worker processes to use when multiprocessing.
N_WORKERS_IN_POOL = os.environ.get('WPHASE_WORKER_COUNT', None)
if N_WORKERS_IN_POOL is not None:
    N_WORKERS_IN_POOL = float(N_WORKERS_IN_POOL)

#: Authority used for wphase.
GA_AUTHORITY = os.environ.get('WPHASE_AUTHORITY', "GA W-phase")

#: The name of the current host.
WPHASE_HOST_NAME = os.environ.get('WPHASE_HOST_NAME', 'localhost')
if not WPHASE_HOST_NAME:
    raise Exception('env var WPHASE_HOST_NAME is defined but is blank')

#: Should profiling information for wphase be produced.
PROFILE_WPHASE = os.environ.get('WPHASE_PROFILE', 'false').lower() in ('1', 'true')

#: Number of traces to put in each wphase results plot
N_TRACES_PER_RESULT_PLOT = os.environ.get('WPHASE_N_TRACES_PER_RESULT_PLOT', 6)

#: Implementation of bandpass filter to use
BANDPASS_IMPLEMENTATION = os.environ.get('WPHASE_BANDPASS_IMPLEMENTATION', 'scipy')

### Model parameters

# Strength of regularization in preliminary magnitude calculation.
# Higher values result in less weight being given to azimuthal variation - in
# the limit at infinity, the azimuths are completely ignored.
# 0.2 seems to work well enough to prevent issues without significantly
# suppressing anisotropy in datasets with good azimuthal coverage.
AMPLITUDE_REGULARIZATION = float(os.environ.get('WPHASE_AMPLITUDE_REGULARIZATION', 0.2))

# Floor to clamp preliminary magnitude to
MINIMUM_PRELIMINARY_MAGNITUDE = float(os.environ.get('WPHASE_MINIMUM_PRELIMINARY_MAGNITUDE', 6.5))

# used to reject a channel if the fitting of the instrument response is
# greater than this.  see "Source inversion of W phase - speeding up
# seismic tsunami warning - Kanamori - 2008" in the doc folder of
# atws_wphase for details on the fitting of the instrument response.
RESPONSE_MISFIT_TOL = 5. # Percent

# time cutoff, seconds/degree from t_p for the Wphase window.
WPHASE_CUTOFF = 15

# traces with peak-to-peak amplitude outside of the range
# [median_rejection_coeff[0] * m, median_rejection_coeff[1] * m], where m
# is the median peak-to-peak amplitude over all stations, will be rejected
MEDIAN_REJECTION_COEFF = [0.1, 3]

# Minimum number of stations required to begin the inversion.
MINIMUM_STATIONS = 4

# Maximum time delay t_d in seconds. We'll search for the optimal t_d up to
# this value.
MAXIMUM_TIME_DELAY = 200.

# Channels where the misfits (100*sqrt(sum(synthetic-observed)^2 /
# sum(synthetic)^2)) are greater than the first element of this sequence  will
# be rejected. If any stations remain, then those with misfits greater than the
# second element will be rejected, and so on.
MISFIT_TOL_SEQUENCE = [300, 200, 100]

# Minimum number of well-fitting channels required for inversion to pass
# quality check
MINIMUM_FITTING_CHANNELS = 10

### Filenames of output products:

#: The name of the JSON file to contain the output from wphase caclcutions.
WPHASE_OUTPUT_FILE_NAME = 'wphase_output.json'

#: Prefix for the names of the images containing the beachball.
WPHASE_BEACHBALL_PREFIX = "beachball"

#: Prefix for the station distribution plot.
WPHASE_STATION_DISTRIBUTION_PREFIX = "station_distribution"

#: Prefix for the grid search plot.
WPHASE_GRID_SEARCH_PREFIX = "grid_search"

#: Prefix for the names of the images containing the traces.
WPHASE_STATION_TRACES_PREFIX = "raw_traces"

#: Prefix for the names of the images containing the wphase traces.
WPHASE_RESULTS_TRACES_PREFIX = "wphase_traces"

#: Prefix for the names of the images containing the wphase misfit info.
WPHASE_MISFIT_PREFIX = "wphase_misfit"

#: Prefix for the names of the images containing the preliminary magnitude fit
WPHASE_PRELIM_FIT_PREFIX = "preliminary_magnitude_fit"


### The remaining settings should probably not be changed - they define the
### structure of the w-phase results JSON.

#: Key for warnings in the result dictionary.
WPHASE_WARNINGS_KEY = 'Warnings'

#: Key for the wphase processing error.
WPHASE_ERROR_KEY = 'Error'

#: Key for the wphase event description.
WPHASE_EVENT_KEY = 'Event'

#: Key containing the stack trace.
WPHASE_ERROR_STACKTRACE_KEY = 'StackTrace'

#: Key for the list of wphase results plots.
RESULTS_PLOTS_KEY = 'WphaseResultPlots'

#: Key for the data source.
WPHASE_DATA_SOURCE_KEY = 'DataSource'

#: Key containing the wpinv profiling output.
WPINV_PROFILE_OUTPUT_KEY = 'WPInvProfile'

#: Key containing the misfits.
MISFIT_KEY = 'Misfits'

#: Key containing the host name.
HOST_NAME_KEY = 'HostName'
