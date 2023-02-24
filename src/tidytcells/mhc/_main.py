'''
Utility functions related to MHCs and MHC genes.
'''


from .._utils.gene_standardisers import HLAStandardiser
from .._utils.standardise_template import standardise_template
from .._utils.warnings import *
import re
from .._resources import HOMOSAPIENS_MHC
from warnings import warn


# --- STATIC RESOURCES ---

CHAIN_ALPHA_RE = re.compile(r'HLA-([ABCEFG]|D[PQR]A)')
CHAIN_BETA_RE = re.compile(r'HLA-D[PQR]B|B2M')
CLASS_1_RE = re.compile(r'HLA-[ABCEFG]|B2M')
CLASS_2_RE = re.compile(r'HLA-D[PQR][AB]')


# --- MAIN FUNCTIONS ---

STANDARDISERS = {
    'homosapiens': HLAStandardiser
}


def standardise(
    gene: str,
    species: str = 'homosapiens',
    precision: str = 'allele'
) -> tuple:
    '''
    Attempt to standardise an MHC gene name to be IMGT-compliant.
    
    An important note here is that this function will verify the validity of an MHC gene up to the level of the protein.
    Any further precise allele designations will not be verified, apart from the requirement that the format (colon-separated numbers) looks valid.
    The reasons for this is firstly because new alleles at that level are added to the IMGT list quite often and so accurate verification is difficult, secondly because people rarely need verification to such a precise level, and finally because such verification costs more computational effort with diminishing returns.

    :param gene:
        Potentially non-standardised MHC gene name.
    :type gene:
        ``str``
    :param species:
        Species to which the MHC gene belongs (see :ref:`supported_species`).
        Defaults to ``'homosapiens'``.
    :type species:
        ``str``
    :param precision:
        The maximum level of precision to standardise to.
        ``'allele'`` standardises to the maximum precision possible.
        ``'protein'`` keeps allele designators up to the level of the protein (first two).
        In this setting, the function returns a tuple instead of a string, where the first element is a string representing the MHC up to the level of the protein,and any further allele designators are separated into a separate string, which is the second element of the tuple (if there are no further designators available, then the second element of the tuple is set to ``None``).
        ``'gene'`` standardises only to the level of the gene.
        Defaults to ``'allele'``.
    :return:
        If the specified ``species`` is supported, and ``gene_name`` could be standardised, then return the standardised gene name.
        If ``species`` is unsupported, then the function does not attempt to standardise, and returns the unaltered ``gene_name`` string.
        Else returns ``None``.
    :rtype:
        ``str``, ``tuple[str]`` (see parameter ``precision``) or ``None``
    '''
    
    # If precision is not either 'allele' or 'gene' raise error
    if not precision in ('allele', 'protein', 'gene'):
        raise ValueError(
            f'precision must be "allele", "protein" or "gene", got {precision}.'
        )
    
    return standardise_template(
        gene=gene,
        gene_type='MHC',
        species=species,
        enforce_functional=True,
        precision=precision,
        standardiser_dict=STANDARDISERS
    )


def get_chain(gene_name: str) -> str:
    '''
    Given a standardised MHC gene name, detect whether it codes for an alpha
    or a beta chain molecule.
    
    :param gene_name:
        Standardised MHC gene name
    :type gene_name:
        ``str``
    :return:
        ``'alpha'`` or ``'beta'`` if ``gene_name`` is recognised and its chain is known, else ``None``.
    :rtype:
        ``str`` or ``None``
    '''

    if type(gene_name) == str:
        # If we don't recognise the gene, return None with warning
        gene_name = gene_name.split('*')[0]

        if not gene_name in (*HOMOSAPIENS_MHC, 'B2M'):
            warn(f'Unrecognised gene {gene_name}. Is this standardised?')
            return None

        if CHAIN_ALPHA_RE.match(gene_name):
            return 'alpha'
        
        if CHAIN_BETA_RE.match(gene_name):
            return 'beta'

        warn(f'Chain for {gene_name} unknown.')
        return None

    raise TypeError(
        f'gene_name must be type str, got {gene_name} ({type(gene_name)}).'
    )


def get_class(gene_name: str) -> int:
    '''
    Given a standardised MHC gene name, detect whether it comprises a class I
    or II MHC receptor complex.
    
    :param gene_name:
        Standardised MHC gene name
    :type gene_name:
        ``str``
    :return:
        ``1`` or ``2`` if ``gene_name`` is recognised and its class is known, else ``None``.
    :rtype:
        ``int`` or ``None``
    '''

    if type(gene_name) == str:
        # If we don't recognise the gene, return None with warning
        gene_name = gene_name.split('*')[0]

        if not gene_name in (*HOMOSAPIENS_MHC, 'B2M'):
            warn(f'Unrecognised gene {gene_name}. Is this standardised?')
            return None

        if CLASS_1_RE.match(gene_name):
            return 1
        
        if CLASS_2_RE.match(gene_name):
            return 2

        warn(f'Class for {gene_name} unknown.')
        return None

    raise TypeError(
        f'gene_name must be type str, got {gene_name} ({type(gene_name)}).'
    )