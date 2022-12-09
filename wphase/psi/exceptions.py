from typing import Optional
from wphase.psi.model import WPhaseResult


class InversionError(Exception):
    """A terminal condition in the W-Phase inversion that is not a malfunction of
    the software itself.

    This is usually raised if there is a problem with input data.

    :param str msg: The reason for termination
    :param WPhaseResult result: Partial result"""
    def __init__(self, *args, result: Optional[WPhaseResult] = None):
        super().__init__(*args)
        self.result = result

class RTdeconvError(Exception):
    """An error encountered during deconvolution of waveforms."""

class UnknownTransferFunction(Exception):
    """An unknown transfer function was encountered in station inventory."""

