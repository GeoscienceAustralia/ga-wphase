# -*- coding: utf-8 -*-

"""
.. module:: settings

Configuration settings for :mod:`wphase`. At present, the settings are declared
directly as python variables, but eventually may be moved to an external file
(e.g. JSON) which persists outside the code base and can be customised for a
specific system.
"""

from __future__ import absolute_import

import os
import sys
import multiprocessing as mp

# Set the number of threads used by OpenBLAS.
os.environ['OPENBLAS_NUM_THREADS'] = '1'

#: Use HDF5 version of Green's functions (the alternative [original] is to use
#: SAC files).
USE_HDF5_GREENS_DB = True

#: Directory containing the green's functions.
GREEN_DIR = os.environ['WPHASE_GREENS_FUNCTIONS']

DEFAULT_OUTPUT_DIR = os.environ.get('WPHASE_OUTPUT_DIR', '/tmp/wphase-output')

#: The number of worker processes to use in :py:class:`multiprocessing.Pool`.
N_WORKERS_IN_POOL = None

#: The name of the JSON file to contain the output from wphase caclcutions.
WPHASE_OUTPUT_FILE_NAME = 'wphase_output.json'

#: Name of the file that contains the exclude list.
WPHASE_EXCLUDE_LIST_FILE_NAME = 'exclude_list.json'

#: wphase intermediate mseed file name.
WPHASE_MSEED_FILE = "Mini-SEED-traces.mseed"

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

#: Authority used for wphase.
GA_AUTHORITY = "GA W-phase"

#: The name of the current host.
WPHASE_HOST_NAME = os.environ.get(
    'WPHASE_HOST_NAME',
    'localhost')
if not WPHASE_HOST_NAME:
    raise Exception('env var WPHASE_HOST_NAME is defined but is blank')

#: Should profiling information for wphase be produced.
PROFILE_WPHASE = False

#: Number of traces to put in each wphase results plot
N_TRACES_PER_RESULT_PLOT = 6

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

#: Implementation of bandpass filter to use
BANDPASS_IMPLEMENTATION = 'scipy'

# Strength of regularization in preliminary magnitude calculation.
# Higher values result in less weight being given to azimuthal variation - in
# the limit at infinity, the azimuths are completely ignored.
# 0.1 seems to work well enough to prevent issues without completely
# suppressing anisotropy.
AMPLITUDE_REGULARIZATION = 0.1

# Floor to clamp preliminary magnitude to
MINIMUM_PRELIMINARY_MAGNITUDE = 6.5
