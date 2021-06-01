"""Simple data classes for representing results internally."""
import numpy as np
from collections import OrderedDict
try:
    from collections.abc import Mapping
except ImportError:
    from collections import Mapping

class Data(Mapping):
    # TODO: switch to dataclasses with real type annotations after ditching
    # python 2
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)

    def __len__(self):
        return len(self.__dict__)
    def __iter__(self):
        return iter(self.__dict__)
    def __getitem__(self, k):
        return self.__dict__[k]

class OL1(Data):
    preliminary_magnitude    = None # type: float
    preliminary_calc_details = None # type: dict

class OL2(OL1):
    moment_tensor            = None # type: np.array
    observed_displacements   = None # type: np.array
    synthetic_displacements  = None # type: np.array
    trace_lengths            = None # type: OrderedDict[str, int]

class OL3(OL2):
    centroid                 = None # type: tuple[float, float, float]
    time_delay               = None # type: float
    grid_search_candidates   = None # type: list
    grid_search_results      = None # type: list[tuple[np.array, float]]
