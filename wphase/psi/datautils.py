from __future__ import absolute_import

from builtins import zip
from builtins import range
import os
from subprocess import PIPE, Popen
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import boxcar,triang, convolve, detrend
from obspy.core import read, Trace, Stream, UTCDateTime

from .bandpass import bandpassfilter
