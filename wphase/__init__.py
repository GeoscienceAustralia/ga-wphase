"""GA W-Phase root module. """

import os
import errno
import logging

from wphase.psi import model

try:
    # to avoid: Exception _tkinter.TclError
    import matplotlib
    matplotlib.use('Agg')
except Exception:
    pass

logger = logging.getLogger(__name__)

from . import settings_schema
settings = settings_schema.WPhaseSettings()

from ._runner_fdsn import runwphase as wphase_runner


def runwphase(
    output_dir = None,
    server = None,
    greens_functions_dir = settings.GREENS_FUNCTIONS,
    n_workers_in_pool = settings.WORKER_COUNT,
    processing_level = 3,
    output_dir_can_exist = False,
    **kwargs) -> model.WPhaseResult:
    """
    Run wphase.

    :param greens_functions_dir: The Green data Directory.
    :param output_dir: Full file path to the output directory. **DO NOT USE
        RELATIVE PATHS**.
    :param n_workers_in_pool: Number of processors to use, (default
        :py:data:`wphase.settings.WORKER_COUNT`) specifies as many as is
        reasonable'.
    :param processing_level: Processing level.
    :param output_dir_can_exist: Can the output directory already exist?
    """

    if server is not None and server.lower() == 'antelope':
        raise Exception('Antelope is no longer supported.')

    # Make the output directory (fails if it already exists).
    if output_dir:
        logger.debug("Creating output directory %s", output_dir)
        try:
            os.makedirs(output_dir)
        except OSError as e:
            if e.errno != errno.EEXIST or not output_dir_can_exist:
                raise

    wphase_results = wphase_runner(
        output_dir,
        server,
        greens_functions_dir,
        n_workers_in_pool,
        processing_level,
        **kwargs)

    wphase_results.HostName = settings.HOST_NAME
    wphase_results.DataSource = server if server else "local files"

    # save the results if output_dir provided
    if output_dir:
        try:
            # TODO: Should this be done in runwphase?
            with open(os.path.join(output_dir,settings.OUTPUT_FILE_NAME), 'w') as of:
                print(wphase_results.json(indent=2), file=of)
        except Exception as e:
            # not sure how we would get here, but we just don't want
            # to stop the rest of processing
            logger.exception("Failed dumping result to JSON.")

    # re-raise any errors from the dark side
    if wphase_results.Error:
        raise Exception(wphase_results.StackTrace)

    return wphase_results

__all__ = ['runwphase']
