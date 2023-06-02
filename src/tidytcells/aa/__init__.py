"""
Functions to clean and standardize amino acid sequence data.
"""


from ._main import standardize


def standardise(*args, **kwargs):
    """
    Alias for :py:func:`tidytcells.aa.standardize`.
    """
    return standardize(*args, **kwargs)
