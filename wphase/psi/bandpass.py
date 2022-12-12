"""Butterworth bandpass filtering"""
from __future__ import absolute_import
from typing import Literal, Union

import numpy as np
from scipy import signal
from wphase import settings

from functools import lru_cache, partial

# The filter design stage is somewhat expensive, so we cache the results
@lru_cache(maxsize=None)
def design_filter(order, nyquist_frequency, low_freq, high_freq):
    corn_freq = [low_freq/nyquist_frequency, high_freq/nyquist_frequency]
    return signal.butter(order, corn_freq, btype='bandpass', output='sos')


def bandpassfilter(
    data: np.ndarray,
    sample_period,
    order,
    low,
    high,
    npass=1,
    impl: Union[Literal["scipy"], Literal["fortran"], None] = None,
    axis: int = -1,
):
    """Apply a bandpass filter along one axis of a numpy array."""
    if not impl:
        impl = settings.BANDPASS_IMPLEMENTATION
    if impl == "scipy":
        sos = design_filter(order, 0.5 / sample_period, low, high)
        return signal.sosfilt(sos, data, axis=axis)
    elif impl == "fortran":
        from . import bpfilter

        if data.ndim > 1:
            return np.apply_along_axis(
                partial(
                    bandpassfilter,
                    sample_period=sample_period,
                    order=order,
                    low=low,
                    high=high,
                    npass=1,
                    impl=impl,
                ),
                axis,
                data,
            )
        return bpfilter.bandpassfilter(
            data, len(data), sample_period, order, npass, low, high
        )
    else:
        raise ValueError("Bandpass implementation %s not known" % impl)
