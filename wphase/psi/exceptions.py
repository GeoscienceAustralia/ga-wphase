class InversionError(Exception):
    """A terminal condition in the W-Phase inversion that is not a malfunction of
    the software itself.

    This is usually raised if there is a problem with input data.

    :param str msg: The reason for termination."""

class RTdeconvError(Exception):
    """An error encountered during deconvolution of waveforms."""
