"""
Functions to clean and standardize MHC gene data.
"""


from ._main import standardize, get_chain, get_class, query
from warnings import warn as _warn


def standardise(*args, **kwargs):
    """
    Alias for :py:func:`tidytcells.mhc.standardize`.
    """
    return standardize(*args, **kwargs)


def classify(*args, **kwargs):
    """
    .. caution:: This will be deprecated soon in favour of :py:func:`tidytcells.mhc.get_class`.

    Alias for :py:func:`tidytcells.mhc.get_class`.
    """
    _warn(
        '"mhc.classify" as an alias will be deprecated in the near future. Please switch to using "mhc.get_class".',
        FutureWarning,
    )
    return get_class(*args, **kwargs)
