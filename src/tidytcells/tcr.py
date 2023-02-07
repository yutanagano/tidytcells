'''
Utility functions related to TCRs and TCR genes.
'''


from .decomposed_gene import _DecomposedGene
import json
from pkg_resources import resource_stream
import re
from typing import Optional
from warnings import warn


# --- STATIC RESOURCES ---


with resource_stream(
    __name__,
    'resources/homosapiens_tcr.json') as s:
    HOMOSAPIENS_TCR = json.load(s)
with resource_stream(
    __name__,
    'resources/homosapiens_tcr_synonyms.json') as s:
    HOMOSAPIENS_TCR_SYNONYMS = json.load(s)

PARSE_RE_HOMOSAPIENS_1 = re.compile(r'^([A-Z0-9\-\.\(\)\/]+)(\*([0-9]+))?$')


# --- HELPER CLASSES ---


class DecomposedTCR(_DecomposedGene):
    def __init__(
        self,
        gene: str,
        allele_designation: Optional[int],
        ref_dict: dict,
        syn_dict: dict
    ) -> None:
        self.gene = gene
        self.allele_designation =\
            f'{allele_designation:02}' if type(allele_designation) == int\
            else None
        self.ref_dict = ref_dict
        self.syn_dict = syn_dict


    @property
    def valid(self) -> bool:
        # Is the gene valid?
        if not self.gene in self.ref_dict:
            return False
        
        # If the gene exists and an allele is specified, check if it exists
        if self.allele_designation:
            return self.allele_designation in self.ref_dict[self.gene]
        
        # If no allele specified, return true as the gene exists
        return True
    

    def compile(self, allele: bool = True) -> str:
        if allele and self.allele_designation:
            return f'{self.gene}*{self.allele_designation}'

        return self.gene
    

    def resolve(self) -> bool:
        # If a synonym, correct to currently approved name
        if self.syn_dict and self.gene in self.syn_dict:
            self.gene = self.syn_dict[self.gene]
        
        return self.valid


# --- HELPER FUNCTIONS ---


def standardise_homosapiens(gene_name: str) -> str:
    # Take note of initial input for reference
    original_input = gene_name

    # Clean whitespace, remove known pollutors
    gene_name = ''.join(gene_name.split())
    gene_name = gene_name.replace('&nbsp;','')

    # Capitalise
    gene_name = gene_name.upper()

    # Parse attempt
    if m := PARSE_RE_HOMOSAPIENS_1.match(gene_name): # ^([A-Z0-9\-\.\(\)\/]+)(\*([0-9]+))?$
        gene = m.group(1)
        allele_designation = None if m.group(3) is None else int(m.group(3))

    # Could not parse
    else:
        warn_failure(original_input, gene_name, 'Homo sapiens')
        return None

    # Build DecomposedTcr object
    decomp_tcr = DecomposedTCR(
        gene=gene,
        allele_designation=allele_designation,
        ref_dict=HOMOSAPIENS_TCR,
        syn_dict=HOMOSAPIENS_TCR_SYNONYMS
    )

    # Try resolving, and return None on failure
    if not decomp_tcr.resolve():
        warn_failure(original_input, decomp_tcr.compile(), 'Homo sapiens')
        return None
    
    return decomp_tcr.compile()


def warn_failure(
    original_input: str,
    attempted_fix: str,
    species: str
):
    warn(
        f'Unrecognised TCR gene name: "{original_input}" '
        f'for species {species}. Attempted fix "{attempted_fix}" did not meet '
        'the standardised format requirements. Ignoring this gene name...'
    )


# --- MAIN FUNCTIONS ---


SUPPORTED_SPECIES = {
    'HomoSapiens': standardise_homosapiens
}


def standardise(gene_name: str, species: str = 'HomoSapiens') -> str:
    '''
    Attempt to standardise a TCR gene name to be IMGT-compliant.

    :param gene_name: Potentially non-standardised TCR gene name.
    :type gene_name: str
    :param species: Species to which the TCR gene belongs (see
        :ref:`supported_species`). Defaults to `"HomoSapiens"`.
    :type species: str
    :return: If the specified ``species`` is supported, and ``gene_name`` could
        be standardised, then return the standardised gene name. If ``species``
        is unsupported, then the function does not attempt to standardise , and
        returns the unaltered ``gene_name`` string. Else return ``None``.
    :rtype: str or None
    
    '''

    # If gene_str is not a string, skip and return None
    if type(gene_name) != str:
        raise TypeError('gene_name must be a str.')

    # If the specified species is supported, attempt standardisation
    if species in SUPPORTED_SPECIES:
        return SUPPORTED_SPECIES[species](gene_name=gene_name)
    
    # Otherwise, don't touch it
    warn(
        f'Unsupported species: "{species}". '
        'Skipping TCR gene standardisation procedure...'
    )
    return gene_name