'Utility functions related to TCRs and TCR genes.'


import json
from pkg_resources import resource_stream
import re
from typing import Union
from warnings import warn


# --- STATIC RESOURCES ---


with resource_stream(
    __name__,
    'resources/tcr_alleles_homosapiens.json') as s:
    _TCR_ALLELES_HOMOSAPIENS = json.load(s)
with resource_stream(
    __name__,
    'resources/tcr_d_designations_homosapiens.json') as s:
    _TCR_D_DESIGNATIONS_HOMOSAPIENS = json.load(s)
with resource_stream(
    __name__,
    'resources/tcr_synonyms_homosapiens.json') as s:
    _TCR_SYNONYMS_HOMOSAPIENS = json.load(s)
PARSE_RE_HOMOSAPIENS_1 = re.compile(
    r'^(TC?R)?([AB][VJ])(\d+)([S-](\d+))?((/|DV|/DV)(\d+))?(/?OR9-2)?([\*\.](.+))?$'
)
PARSE_RE_HOMOSAPIENS_2 = re.compile(
    r'^(\d+)([S-](\d+))?(/?DV(\d+))(\*(\d+))?$'
)


# --- HELPER CLASSES ---


class DecomposedTcr:
    def __init__(
        self,
        base: Union[str, None],
        num1: Union[str, None],
        num2: Union[str, None],
        p: bool,
        or92: bool,
        d_designation: Union[str, None],
        allele_num: Union[str, None]
    ) -> None:
        self.base = base
        self.num1 = num1
        self.num2 = num2
        self.p = p
        self.or92 = or92
        self.d_designation = d_designation
        self.allele_num = allele_num


    @property
    def valid(self) -> bool:
        gene_str = self.compile(allele=False)

        # Is the gene valid?
        if not gene_str in _TCR_ALLELES_HOMOSAPIENS:
            return False
        
        # If the gene exists and an allele is specified, check if it exists
        if self.allele_num:
            return self.allele_num in _TCR_ALLELES_HOMOSAPIENS[gene_str]
        
        # If no allele specified, return true as the gene exists
        return True
    

    def compile(self, allele: bool = True) -> str:
        compiled = 'TR' + self.base + self.num1

        if self.num2:
            compiled = compiled + '-' + self.num2

        if self.p:
            compiled = compiled + 'P'

        if self.d_designation:
            compiled = compiled + 'DV' + self.d_designation
        
        if self.or92:
            compiled = compiled + 'OR9-2'
        
        if allele and self.allele_num:
            compiled = compiled + '*' + self.allele_num

        return compiled


    def resolve(self) -> None:
        # Remove any leading zeros in gene numbers
        self.num1 = None if self.num1 is None else self.num1.lstrip('0')
        self.num2 = None if self.num2 is None else self.num2.lstrip('0')

        # Resolve TRD designation if appropriate
        if self.d_designation is None:
            try:
                self.d_designation = _TCR_D_DESIGNATIONS_HOMOSAPIENS[
                    self.compile(allele=False)
                ]

            except(KeyError):
                pass

        # Add leading zero to allele number if appropriate
        self._resolve_allele_num()

        # If now a valid TCR, return True
        if self.valid:
            return True

        # If still invalid, try nullifying num2 if appropriate
        if self.num2 == '1':
            self.num2 = None

            if self.valid:
                return True
            
            self.num2 = '1'
        
        return False


    def _resolve_allele_num(self) -> None:
        '''
        Either clean up allele_num if it contains something that looks like a
        real number, otherwise set it to None.
        '''
        if self.allele_num is None:
            return

        if not self.allele_num.isdigit():
            self.allele_num = None
            return

        if len(self.allele_num) == 1:
            self.allele_num = '0' + self.allele_num


# --- HELPER FUNCTIONS ---


def _standardise_homosapiens(gene_name: str) -> str:
    # Take note of initial input for reference
    original_input = gene_name

    # Clean whitespace, remove known pollutors
    gene_name = ''.join(gene_name.split())
    gene_name = gene_name.replace('&nbsp;','')

    # Parse attempt 1
    if m := PARSE_RE_HOMOSAPIENS_1.match(gene_name):
        gene_prefix = m.group(1)
        base = m.group(2)
        num1 = m.group(3)
        num2 = m.group(5)
        d_designation = m.group(8)
        or92 = True if m.group(9) else False
        allele_num = m.group(11)

        # If input looks like a deprecated name, see if it can be resolved to
        # its modern couterpart via our synonym lookup table
        if gene_prefix == 'TCR':
            reconstructed = 'TCR' + base + num1
            if num2:
                reconstructed = reconstructed + 'S' + num2

            try:
                translation = _TCR_SYNONYMS_HOMOSAPIENS[reconstructed]
                base = translation['base']
                num1 = translation['num1']
                num2 = translation['num2']
                or92 = translation['OR9-2']

            except(KeyError):
                pass

    # Parse attempt 2
    elif m := PARSE_RE_HOMOSAPIENS_2.match(gene_name):
        gene_prefix = None
        base = 'AV'
        num1 = m.group(1)
        num2 = m.group(3)
        d_designation = m.group(5)
        or92 = False
        allele_num = m.group(7)

    # Could not parse
    else:
        _warn_failure(original_input, gene_name, 'Homo sapiens')
        return None

    # Build DecomposedTcr object
    decomp_tcr = DecomposedTcr(
        base=base,
        num1=num1,
        num2=num2,
        p=False,
        or92=or92,
        d_designation=d_designation,
        allele_num=allele_num
    )

    # Try resolving, and return None on failure
    if not decomp_tcr.resolve():
        _warn_failure(original_input, decomp_tcr.compile(), 'Homo sapiens')
        return None
    
    return decomp_tcr.compile()


def _warn_failure(
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
    'HomoSapiens': _standardise_homosapiens
}


def standardise(gene_name: str, species: str) -> str:
    # If gene_str is not a string, skip and return None
    if type(gene_name) != str:
        _warn_failure(gene_name, gene_name, species)
        return None

    # If the specified species is supported, attempt standardisation
    if species in SUPPORTED_SPECIES:
        return SUPPORTED_SPECIES[species](gene_name=gene_name)
    
    # Otherwise, don't touch it
    warn(
        f'Unsupported species: "{species}". '
        'Skipping TCR gene standardisation procedure...'
    )
    return gene_name