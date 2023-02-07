'''
Utility functions related to MHCs and MHC genes.
'''


from .decomposed_gene import _DecomposedGene
from itertools import product
import json
from pkg_resources import resource_stream
import re
from typing import Tuple, Optional
from warnings import warn


# --- STATIC RESOURCES ---


with resource_stream(__name__, 'resources/homosapiens_mhc.json') as s:
    HOMOSAPIENS_MHC = json.load(s)
with resource_stream(__name__, 'resources/homosapiens_mhc_synonyms.json') as s:
    HOMOSAPIENS_MHC_SYNONYMS = json.load(s)

PARSE_RE = re.compile(r'^([A-Z0-9\-]+)(\*([\d:]+G?P?)[LSCAQN]?)?$')


# --- HELPER CLASSES ---


class DecomposedMHC(_DecomposedGene):
    def __init__(
        self,
        gene: str,
        allele_designation: Optional[list[str]],
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
    

    def compile(self, allele: bool = True) -> str:
        if allele and self.allele_designation:
            return f'{self.gene}*{":".join(self.allele_designation)}'

        return self.gene
    

    def resolve(self) -> bool:
        # If a synonym, correct to currently approved name
        if self.syn_dict and self.gene in self.syn_dict:
            self.gene = self.syn_dict[self.gene]
        
        return self.valid


# --- HELPER FUNCTIONS ---


def warn_failure(
    original_input: str,
    attempted_fix: str,
    species: str
) -> None:
    warn(
        f'Unrecognised MHC gene name: "{original_input}" '
        f'for species {species}. Attempted fix "{attempted_fix}" did not meet '
        'the standardised format requirements. Ignoring this gene name...'
    )


# --- MAIN FUNCTIONS ---


SUPPORTED_SPECIES = {
    'HomoSapiens': (HOMOSAPIENS_MHC, HOMOSAPIENS_MHC_SYNONYMS)
}


def standardise(gene_name: str, species: str = 'HomoSapiens') -> tuple:
    '''
    Attempt to standardise an MHC gene name to be IMGT-compliant.

    :param gene_name: Potentially non-standardised MHC gene name.
    :type gene_name: str
    :param species: Species to which the MHC gene belongs (see
        :ref:`supported_species`). Defaults to `"HomoSapiens"`.
    :type species: str
    :return: If the specified ``species`` is supported, and ``gene_name`` could
        be standardised, then return a tuple containing the standardised gene
        name decomposed into two parts: 1) the name of the gene specific to the
        level of the protein, and 2) (if any) further valid specifier fields
        (e.g. ``('HLA-A*01:01', ':01:01')``, see :ref:`example_usage`). If
        ``species`` is unsupported, then the function does not attempt to
        standardise, and returns a tuple with the unaltered ``gene_name`` for
        the first element, and ``None`` for the second element. Else return the
        tuple ``(None, None)``.
    :rtype: tuple[str or None] or None

    '''

    # If gene_str is not a string, skip and return None.
    if type(gene_name) != str:
        raise TypeError('gene_name must be a str.')

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
    if m := PARSE_RE.match(gene_name): # ^^([A-Z0-9\-]+)(\*([\d:]+G?P?)[LSCAQN]?)?$
        gene = m.group(1)
        allele_designation = None if m.group(3) is None else m.group(3).split(':')

    # Could not parse
    else:
        warn_failure(original_input, gene_name, species)
        return None

    # Build DecomposedMHC object
    decomp_mhc = DecomposedMHC(
        gene=gene,
        allele_designation=allele_designation,
        ref_dict=ref_dict,
        syn_dict=syn_dict
    )

    # Try resolving, and return None on failure
    if not decomp_mhc.resolve():
        warn_failure(original_input, decomp_mhc.compile(), species)
        return None
    
    return decomp_mhc.compile()