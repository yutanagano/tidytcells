"""
Functions to clean and standardise junction (CDR3) data.
"""


from ._main import standardise


def standardize(*args, **kwargs):
    """
    Alias for :py:func:`tidytcells.junction.standardise`.
    """
    return standardise(*args, **kwargs)
