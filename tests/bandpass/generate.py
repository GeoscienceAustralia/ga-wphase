#!/usr/bin/python
# Generates test datasets for bandpass filtering.
# This was originally used to generate the test files that are now checked in
# to the repository, using impl="fortran".  Thus test_bandpass now verifies
# that the new impl="scipy" method produces the same results.

import os
import sys
import numpy as np
import json
from wphase.psi.bandpass import bandpassfilter

PARAM_SETS = [
    (1.0, 1.0,  4, 0.002, 0.008),
    (0.5, 10.0, 4, 0.002, 0.008),
    (1.0, 0.1,  8, 0.002, 0.008),
    (1.0, 10.0, 8, 0.002, 0.02),
    (0.1, 1.0,  8, 0.02,  0.1),
    (5.0, 1.0,  4, 0.001, 0.003),
]

def generate_test_data(impl):
    for i, params in enumerate(PARAM_SETS):
        generate(impl, str(i), *params)

def generate(impl, dirname, delta, amp, order, lowfreq, hifreq):
    data = np.random.normal(0, amp, 5000)
    filtered = bandpassfilter(data, delta, order, lowfreq, hifreq, impl=impl)
    try:
        os.mkdir(dirname)
    except Exception:
        pass
    np.save(dirname + '/before.npy', data)
    np.save(dirname + '/after.npy', filtered)
    with open(dirname + '/params.json', 'wt') as fh:
        json.dump(dict(
            delta=delta,
            order=order,
            lowfreq=lowfreq,
            hifreq=hifreq,
            amp=amp), fh)

if __name__ == '__main__':
    impl = sys.argv[1]
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    generate_test_data(impl)
