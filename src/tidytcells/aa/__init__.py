"""
Functions to clean and standardise amino acid sequence data.
"""


from ._main import standardise


def standardize(*args, **kwargs):
    """
    Alias for :py:func:`tidytcells.aa.standardise`.
    """
    return standardise(*args, **kwargs)
