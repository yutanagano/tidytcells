'''
Utility functions related to TCRs and TCR genes.
'''


from .._utils.gene_standardisers import (
    HomoSapiensTCRStandardiser,
    MusMusculusTCRStandardiser
)
from .._utils.warnings import *


STANDARDISERS = {
    'homosapiens': HomoSapiensTCRStandardiser,
    'musmusculus': MusMusculusTCRStandardiser
}


def standardise(
    gene: str,
    species: str = 'homosapiens',
    enforce_functional: bool = False,
    precision: str = 'allele'
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
    :return:
        If the specified ``species`` is supported, and ``gene_name`` could be standardised, then return the standardised gene name.
        If ``species`` is unsupported, then the function does not attempt to standardise , and returns the unaltered ``gene_name`` string.
        Else returns ``None``.
    :rtype:
        ``str`` or ``None``
    '''

    # If gene_str is not a string, raise error
    if type(gene) != str:
        raise TypeError(
            f'gene_name must be type str, got '
            f'{gene} ({type(gene)}).'
        )
    
    # If precision is not either 'allele' or 'gene' raise error
    if not precision in ('allele', 'gene'):
        raise ValueError(
            f'precision must be either "allele" or "gene", got {precision}.'
        )

    # If the specified species is not supported, no-op (with warning)
    if not species in STANDARDISERS:
        warn_unsupported_species(species, 'TCR')
        return gene

    # Take note of initial input for reference
    original_input = gene

    # Clean whitespace, remove known pollutors, uppercase
    gene = ''.join(gene.split())
    gene = gene.replace('&nbsp;','')
    gene = gene.upper()

    # Standardisation attempt
    standardised = STANDARDISERS[species](gene)

    if not standardised.valid(enforce_functional): # Standaridsation failure
        warn_failure(original_input, standardised.compile(precision), species)
        return None
    
    return standardised.compile(precision)