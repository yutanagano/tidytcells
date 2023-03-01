'''
Utility functions related to TCRs and TCR genes.
'''


from typing import Optional
from .._utils.gene_standardisers import (
    HomoSapiensTCRStandardiser,
    MusMusculusTCRStandardiser
)
from .._utils.standardise_template import standardise_template
from .._utils.warnings import *


STANDARDISERS = {
    'homosapiens': HomoSapiensTCRStandardiser,
    'musmusculus': MusMusculusTCRStandardiser
}


def standardise(
    gene: Optional[str] = None,
    species: str = 'homosapiens',
    enforce_functional: bool = False,
    precision: str = 'allele',

    gene_name: Optional[str] = None
) -> str:
    '''
    Attempt to standardise a TCR gene name to be IMGT-compliant.

    :param gene:
        Potentially non-standardised TCR gene name.
    :type gene:
        ``str``
    :param species:
        Species to which the TCR gene belongs (see :ref:`supported_species`).
        Defaults to ``'homosapiens'``.
    :type species:
        ``str``
    :param enforce_functional:
        If ``True``, disallows TCR genes that are recognised by IMGT but are marked as non-functional (ORF or pseudogene).
        Defaults to ``False``.
    :type enforce_functional:
        ``bool``
    :param precision:
        The maximum level of precision to standardise to.
        ``'allele'`` standardises to the maximum precision possible.
        ``'gene'`` standardises only to the level of the gene.
        Defaults to ``'allele'``.
    :type precision:
        ``str``

    :param gene_name:
        Alias for the parameter ``gene``.
    :type gene_name:
        ``str``

    :return:
        If the specified ``species`` is supported, and ``gene`` could be standardised, then return the standardised gene name.
        If ``species`` is unsupported, then the function does not attempt to standardise , and returns the unaltered ``gene`` string.
        Else returns ``None``.
    :rtype:
        ``str`` or ``None``
    '''
    # Alias resolution
    if gene is None:
        gene = gene_name
    
    # If precision is not either 'allele' or 'gene' raise error
    if not precision in ('allele', 'gene'):
        raise ValueError(
            f'precision must be either "allele" or "gene", got {precision}.'
        )

    return standardise_template(
        gene=gene,
        gene_type='TCR',
        species=species,
        enforce_functional=enforce_functional,
        precision=precision,
        standardiser_dict=STANDARDISERS
    )