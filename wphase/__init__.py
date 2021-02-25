"""
wphase calculations.
"""

try:
    # to avoid: Exception _tkinter.TclError
    import matplotlib
    matplotlib.use('Agg')
except Exception:
    pass

import os
import json

import wphase.settings
from wphase._runner_fdsn import runwphase as wphase_runner

def runwphase(
    output_dir,
    server,
    greens_functions_dir = settings.GREEN_DIR,
    n_workers_in_pool = settings.N_WORKERS_IN_POOL,
    processing_level = 3,
    stations_to_exclude = None,
    output_dir_can_exist = False,
    **kwargs):
    """
    Run wphase.

    :param greens_functions_dir: The Green data Directory.
    :param output_dir: Full file path to the output directory. **DO NOT USE
        RELATIVE PATHS**.
    :param n_workers_in_pool: Number of processors to use, (default
        :py:data:`wphase.settings.N_WORKERS_IN_POOL`) specifies as many as is
        reasonable'.
    :param processing_level: Processing level.
    :param stations_to_exclude: List of station identifiers to exclude.
    :param output_dir_can_exist: Can the output directory already exist?
    """

    if server.lower() == 'antelope':
        raise Exception('Antelope is no longer supported.')

    # Make the output directory (fails if it already exists).
    if output_dir_can_exist:
        try: os.makedirs(output_dir)
        except OSError: pass
    else:
        os.makedirs(output_dir)

    wphase_results = wphase_runner(
        output_dir,
        server,
        greens_functions_dir,
        n_workers_in_pool,
        processing_level,
        stations_to_exclude,
        output_dir_can_exist,
        **kwargs)

    wphase_results[settings.HOST_NAME_KEY] = settings.WPHASE_HOST_NAME
    wphase_results[settings.WPHASE_DATA_SOURCE_KEY] = server

    # save the results
    with open(os.path.join(output_dir, settings.WPHASE_OUTPUT_FILE_NAME), 'w') as out:
        json.dump(wphase_results, out)

    # re-raise any errors from the dark side
    if settings.WPHASE_ERROR_KEY in wphase_results:
        raise Exception(wphase_results[settings.WPHASE_ERROR_STACKTRACE_KEY])

    return wphase_results
