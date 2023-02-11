'''
Utility functions related to MHCs and MHC genes.
'''


from .._decomposed_gene import _DecomposedGene
import json
from pkg_resources import resource_stream as rs
import re
from typing import List, Optional
from warnings import warn


# --- STATIC RESOURCES ---


with rs('tidytcells', 'resources/homosapiens_mhc.json') as s:
    HOMOSAPIENS_MHC = json.load(s)
with rs('tidytcells', 'resources/homosapiens_mhc_synonyms.json') as s:
    HOMOSAPIENS_MHC_SYNONYMS = json.load(s)

PARSE_RE = re.compile(r'^([A-Z0-9\-\.\:\/]+)(\*([\d:]+G?P?)[LSCAQN]?)?')

FIX_RE_1 = re.compile(r'^([A-Z\-\.\:\/]+)(\d+)$')

CHAIN_ALPHA_RE = re.compile(r'HLA-([ABCEFG]|D[PQR]A)')
CHAIN_BETA_RE = re.compile(r'HLA-D[PQR]B|B2M')
CLASS_1_RE = re.compile(r'HLA-[ABCEFG]|B2M')
CLASS_2_RE = re.compile(r'HLA-D[PQR][AB]')


# --- HELPER CLASSES ---


class DecomposedHLA(_DecomposedGene):
    def __init__(
        self,
        gene: str,
        allele_designation: Optional[List[str]],
        precision: str,
        ref_dict: dict,
        syn_dict: dict
    ) -> None:
        self.gene = gene
        self.allele_designation =\
            [] if allele_designation is None\
            else allele_designation
        if not allele_designation is None:
            self.is_g = allele_designation[-1].endswith('G')
            self.is_p = allele_designation[-1].endswith('P')
        else:
            self.is_g = self.is_p = False
        self.precision = precision

        self.ref_dict = ref_dict
        self.syn_dict = syn_dict


    @property
    def valid(self) -> bool:
        # Is the gene valid?
        if not self.gene in self.ref_dict:
            return False

        # If the gene exists and there are allele designators, walk down
        # allele tree up to the level of the protein (or the level of the G or
        # P group when appropriate) to see if the designator field values are
        # valid
        allele_designation = self.allele_designation.copy()
        if not (self.is_g or self.is_p):
            allele_designation = allele_designation[:2]
        current_root = self.ref_dict[self.gene]

        while len(allele_designation) > 0:
            try:
                current_root = current_root[allele_designation.pop(0)]
            except(KeyError):
                # We exit the tree at some point, so designator fields are
                # invalid
                return False

        # If there are designator fields past the protein level, just make sure
        # they look like legitimate designator field values
        if not (self.is_g or self.is_p) and len(self.allele_designation) > 2:
            further_designators = self.allele_designation[2:]

            if len(further_designators) > 2:
                return False

            for field in further_designators:
                if not field.isdigit():
                    return False
                
                if len(field) < 2:
                    return False
        
        # Valid gene with valid specifier fields, so valid!
        return True
    

    def compile(self) -> str:
        if self.allele_designation:
            if self.precision == 'allele':
                return f'{self.gene}*{":".join(self.allele_designation)}'
            
            if self.precision == 'protein':
                prot_designation = self.allele_designation[:2]
                further_designation =\
                    None if len(self.allele_designation) <= 2\
                    else ':'+':'.join(self.allele_designation[2:])
                
                return (
                    f'{self.gene}*{":".join(prot_designation)}',
                    further_designation
                )

        return self.gene
    

    def resolve(self) -> bool:
        # If a synonym, correct to currently approved name
        if self.syn_dict and self.gene in self.syn_dict:
            self.gene = self.syn_dict[self.gene]
        
        if not self.gene.startswith('HLA-'):
            self.gene = 'HLA-' + self.gene
        
        if self.valid:
            return True
        
        m = FIX_RE_1.match(self.gene)
        if not self.allele_designation and m:
            self.gene = m.group(1)
            self.allele_designation = [f'{int(m.group(2)):02}']

        return self.valid


# --- HELPER FUNCTIONS ---


def warn_failure(
    original_input: str,
    attempted_fix: str,
    species: str
) -> None:
    warn(
        f'Failed to standardise: "{original_input}" for species {species}. '
        f'Attempted fix "{attempted_fix}" did not meet the standardised '
        'format requirements. Ignoring this gene name...'
    )


# --- MAIN FUNCTIONS ---


SUPPORTED_SPECIES = {
    'HomoSapiens': (HOMOSAPIENS_MHC, HOMOSAPIENS_MHC_SYNONYMS)
}


def standardise(
    gene_name: str,
    species: str = 'HomoSapiens',
    precision: str = 'allele'
) -> tuple:
    '''
    Attempt to standardise an MHC gene name to be IMGT-compliant.

    :param gene_name:
        Potentially non-standardised MHC gene name.
    :type gene_name:
        ``str``
    :param species:
        Species to which the MHC gene belongs (see :ref:`supported_species`).
        Defaults to ``'HomoSapiens'``.
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

    # If gene_str is not a string, skip and return None.
    if type(gene_name) != str:
        raise TypeError(
            f'gene_name must be type str, got '
            f'{gene_name} ({type(gene_name)}).'
        )
    
    # If precision is not either 'allele' or 'gene' raise error
    if not precision in ('allele', 'protein', 'gene'):
        raise ValueError(
            f'precision must be "allele", "protein" or "gene", got {precision}.'
        )

    # If the specified species is not supported, no-op (with warning)
    if not species in SUPPORTED_SPECIES:
        # Otherwise, don't touch it
        warn(
            f'Unsupported species: "{species}". '
            'Skipping MHC gene standardisation procedure...'
        )
        return gene_name
    
    ref_dict, syn_dict = SUPPORTED_SPECIES[species]

    # Take note of initial input for reference
    original_input = gene_name

    # Clean whitespace, remove known pollutors
    gene_name = ''.join(gene_name.split())
    gene_name = gene_name.replace('&nbsp;','')

    # Capitalise
    gene_name = gene_name.upper()

    # Return B2M as is
    if gene_name == 'B2M':
        return 'B2M'

    # Parse attempt
    m = PARSE_RE.match(gene_name) # ^([A-Z0-9\-\.\:\/]+)(\*([\d:]+G?P?)[LSCAQN]?)?$
    if m:
        gene = m.group(1)
        allele_designation =\
            None if m.group(3) is None else m.group(3).split(':')

    # Could not parse
    else:
        warn_failure(original_input, gene_name, species)
        return None

    # Build DecomposedMHC object
    decomp_mhc = DecomposedHLA(
        gene=gene,
        allele_designation=allele_designation,
        precision=precision,
        ref_dict=ref_dict,
        syn_dict=syn_dict
    )

    # Try resolving, and return None on failure
    if not decomp_mhc.resolve():
        warn_failure(original_input, decomp_mhc.compile(), species)
        return None
    
    return decomp_mhc.compile()


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