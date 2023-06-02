"""
Functions to clean and standardize MHC gene data.
"""


from ._main import standardize, get_chain, get_class, query


def standardise(*args, **kwargs):
    """
    Alias for :py:func:`tidytcells.mhc.standardize`.
    """
    return standardize(*args, **kwargs)


def classify(*args, **kwargs):
    """
    Alias for :py:func:`tidytcells.mhc.get_class`.
    """
    return get_class(*args, **kwargs)
