"""Greens function data source."""
import logging
import os.path

from future.utils import raise_from

import h5py
import numpy as np
import obspy

logger = logging.getLogger(__name__)

try:
    from functools import lru_cache
except ImportError:
    from backports.functools_lru_cache import lru_cache

# We follow Kanamori & Rivera in using spherical coordinates r, theta, phi
# where r is vertical, theta is latitude and phi is longitude.

# Moment tensors are symmetric 3x3 matrices, so we store them as 6x1 vectors
# with the following order convention (R=r, P=phi, T=theta):
ALL_MT_COMPONENTS = ["RR", "PP", "TT", "TP", "RT", "RP"]

# Thanks to symmetry, we can reduce any situation to one where the station is
# due north of the source, which then implies that certain combinations of
# waveform component and moment tensor component are always zero. This data
# structure stores the non-zero combinations.
NONZERO_COMPONENTS = dict(
    Z=["PP", "RR", "RT", "TT"],
    L=["PP", "RR", "RT", "TT"],
    T=["TP", "RP"],
)

# East/North is achieved by rotating Transverse/Longitudinal
SYMMETRIZED_COMPONENT = {"Z": "Z", "N": "L", "E": "T"}

class GreensFunctions(object):
    """Class that provides convenient access to the Greens Function database
    used by Kanamori & Rivera 2008.

    This can be stored either in a directory structure or a HDF5 file."""
    dist_scale = 10.
    dist_offset = 0.5
    dist_max = 90.

    # greens functions are not stored in standard units - we multiply by this
    # factor to convert to Newton-metres:
    gf_unit_Nm = 10.**-21

    def __init__(self, path):
        """:param str path: Path to the GF directory *or* HDF5 file."""
        self.path = path
        self.is_hdf5 = not os.path.isdir(path)
        if self.is_hdf5:
            with self._hdf_open() as f:
                self._hdirs = list(f)
        else:
            self._hdirs = [f for f in os.listdir(path)
                           if f.endswith('.5')]
        self.depths = np.array([float(d[1:]) for d in self._hdirs], dtype=float)

    def _hdf_open(self):
        return h5py.File(self.path, 'r')

    @lru_cache(maxsize=None)
    def _closest_depth(self, depth):
        idx = np.argmin(np.abs(self.depths - depth))
        return self._hdirs[idx] + '/'

    def _get_array(self, path):
        if self.is_hdf5:
            with self._hdf_open() as f:
                return f[path][...]
        else:
            return self._read_obspy_trace(path)[0].data

    def _read_obspy_trace(self, path):
        return obspy.read(os.path.join(self.path, path))

    def __getitem__(self, tup):
        hdir, dist, mt_component, wf_component = tup
        if wf_component not in SYMMETRIZED_COMPONENT.values():
            raise ValueError("Bad component %s" % wf_component)
        try:
            return self._get_array(
                "%s%s/GF.%s.SY.LH%s.SAC" %
                (hdir, mt_component, dist, wf_component)
            ) * self.gf_unit_Nm
        except Exception as e:
            raise_from(KeyError(tup), e)

    def select(self, wf_component, dist, depth):
        """Retrieves the raw Greens functions for the given distance, depth and
        component.

        :param float dist:
            The epicentral distance in degrees.
        :param float depth:
            The depth.
        :param str wf_component:
            Waveform component in reduced coordinates: Either Z, T or L.
        :return:
            6xN array containing the six greens functions corresponding to the
            requested waveform component and each of the six moment tensor components.
            (Some of these will be zero due to symmetry.)
        :rtype:
            np.array
        """
        hdir = self._closest_depth(depth)

        dist = min(dist, self.dist_max)
        dist_form = str(int(dist*self.dist_scale+self.dist_offset)).zfill(4)

        mt_components = NONZERO_COMPONENTS[wf_component]

        rows = [np.asarray(self[hdir, dist_form, mt_component, wf_component], dtype=np.float64)
                if mt_component in mt_components
                else None
                for mt_component in ALL_MT_COMPONENTS]
        shape = next(r.shape for r in rows if r is not None)
        return np.array([r if r is not None
                         else np.zeros(shape, dtype=np.float64)
                         for r in rows])

    def select_rotated(self, wf_component, dist, depth, azimuth):
        """Retrieves the transformed Greens functions for the given azimuth,
        distance, depth and component.

        :param float dist:
            The epicentral distance in degrees.
        :param float depth:
            The depth.
        :param str wf_component:
            Waveform component in world coordinates: Either Z, E or N.
        :param float azi:
            Azimuth in radians.
        :return:
            6xN array containing the six greens functions corresponding to the
            requested waveform component and each of the six moment tensor components.
            (Some of these will be zero due to symmetry.)
        :rtype:
            np.array
        """
        if wf_component == 'Z':
            # vertical waveform component, so no waveform rotation necessary
            unrotated = self.select('Z', dist, depth)
            return GFMTtransform(azimuth, wf_component).dot(unrotated)
        else:
            L = self.select('L', dist, depth)
            T = self.select('T', dist, depth)

            # First transform moment tensor components:
            Ld = GFMTtransform(azimuth, 'N').dot(L)
            Td = GFMTtransform(azimuth, 'E').dot(T)

            # Then waveform components:
            if wf_component in ('N', '1'):
                return Ld * np.cos(azimuth) + Td * np.sin(azimuth)
            elif wf_component in ('E', '2'):
                return Ld * np.sin(azimuth) - Td * np.cos(azimuth)
            else:
                raise ValueError("Unknown waveform component %s" % wf_component)


    @property
    @lru_cache()
    def delta(self):
        """Time resolution of greens functions database"""
        if self.is_hdf5:
            with self._hdf_open() as f:
                return f.attrs["dt"]
        else:
            return self._read_obspy_trace("H003.5/PP/GF.0001.SY.LHZ.SAC")[0].stats.delta

@lru_cache(maxsize=None)
def GFMTtransform(azimuth, wf_component):
    """Build the transformation matrix corresponding to a rotation by the given
    azimuth that can be applied to moment tensor greens functions for the given
    waveform component."""
    s = np.sin(azimuth)
    c = np.cos(azimuth)
    # ["RR", "PP", "TT", "TP", "RT", "RP"]
    if wf_component in ('Z', 'N'):
        return np.array([
            [1, 0,       0,      0, 0, 0],
            [0, c*c,     s*s,    0, 0, 0],
            [0, s*s,     c*c,    0, 0, 0],
            [0, -2.*s*c, 2.*s*c, 0, 0, 0],
            [0, 0,       0,      0, c, 0],
            [0, 0,       0,      0, s, 0],
        ])
    elif wf_component == 'E':
        s2 = np.sin(2.*azimuth)
        c2 = np.cos(2.*azimuth)
        return np.array([
            [0, 0, 0, 0,       0, 0],
            [0, 0, 0, 0.5*s2,  0, 0],
            [0, 0, 0, -0.5*s2, 0, 0],
            [0, 0, 0, c2,      0, 0],
            [0, 0, 0, 0,       0, -s],
            [0, 0, 0, 0,       0, c],
        ])


# SelectN => PP, RR, RT, TT
# SelectE => RP, TP
# SelectZ => PP, RR, RT, TT
