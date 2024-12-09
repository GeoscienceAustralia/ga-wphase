"""Butterworth bandpass filtering"""

from scipy import signal
from wphase import settings

from functools import lru_cache

# The filter design stage is somewhat expensive, so we cache the results
@lru_cache(maxsize=None)
def design_filter(order, nyquist_frequency, low_freq, high_freq):
    corn_freq = [low_freq/nyquist_frequency, high_freq/nyquist_frequency]
    return signal.butter(order, corn_freq, btype='bandpass', output='sos')

def bandpassfilter(data, sample_period, order, low, high, npass=1, impl=None):
    if not impl:
        impl = settings.BANDPASS_IMPLEMENTATION
    if impl == 'scipy':
        sos = design_filter(order, 0.5 / sample_period, low, high)
        return signal.sosfilt(sos, data)
    elif impl == 'fortran':
        from . import bpfilter
        return bpfilter.bandpassfilter(data, len(data), sample_period,
                                       order, npass, low, high)
    else:
        raise ValueError("Bandpass implementation %s not known" % impl)
