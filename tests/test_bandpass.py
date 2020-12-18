import json
import os
from glob import glob

import numpy as np
import pytest
from numpy.linalg import norm
from wphase.psi.bandpass import bandpassfilter

# Read test datasets from disk.
# These were generated using bandpass/generate.py and the old SAC filter.
TESTDIR = os.path.join(os.path.dirname(__file__), 'bandpass')
CASEDIRS = glob(TESTDIR + '/*/')
def read_case(casedir):
    return (
        json.load(open(casedir + '/params.json')),
        np.load(casedir + '/before.npy'),
        np.load(casedir + '/after.npy'),
    )

TOLERANCE_FACTOR = 1e-4

@pytest.mark.parametrize('params,before,after', map(read_case, CASEDIRS))
def test_scipy_produces_same_results(params, before, after):
    """We replaced a fortran implementation of a second-order-section band-pass
    filter with scipy's sos_filt. This test is to ensure that the two
    implementations produce the same results."""
    delta = params['delta']
    order = params['order']
    lowfreq = params['lowfreq']
    hifreq = params['hifreq']
    filtered = bandpassfilter(before, delta, order, lowfreq, hifreq)
    assert norm(filtered - after) < TOLERANCE_FACTOR * norm(before - after)
