"""
Functions to clean and standardise TCR gene data.
"""


from ._main import standardize, query, get_aa_sequence


def standardise(*args, **kwargs):
    """
    Alias for :py:func:`tidytcells.tcr.standardize`.
    """
    return standardize(*args, **kwargs)
