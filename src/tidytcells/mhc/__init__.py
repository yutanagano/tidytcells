'''
Functions to clean and standardise MHC gene data.
'''


from ._main import (
    standardise,
    get_chain,
    get_class
)


def standardize(*args, **kwargs):
    '''
    Alias for :py:func:`tidytcells.mhc.standardise`.
    '''
    return standardise(*args, **kwargs)