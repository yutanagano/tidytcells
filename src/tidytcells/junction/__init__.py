"""
Functions to clean and standardize junction (CDR3) data.
"""


from ._main import standardize


def standardise(*args, **kwargs):
    """
    Alias for :py:func:`tidytcells.junction.standardize`.
    """
    return standardize(*args, **kwargs)
