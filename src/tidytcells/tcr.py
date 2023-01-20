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
    r'^(\d+)([S-](\d+))?(/?DV(\d+))([\*\.](\d+))?$'
)


with resource_stream(
    __name__,
    'resources/tcr_alleles_musmusculus.json') as s:
    _TCR_ALLELES_MUSMUSCULUS = json.load(s)
with resource_stream(
    __name__,
    'resources/tcr_d_designations_musmusculus.json') as s:
    _TCR_D_DESIGNATIONS_MUSMUSCULUS = json.load(s)

PARSE_RE_MUSMUSCULUS_1 = re.compile(
    r'^TR([AB][VJ])(\d+)(D)?([S-](\d+))?((/|DV|/DV)(\d+)(D)?(-(\d))?)?([\*\.](.+))?$'
)
PARSE_RE_MUSMUSCULUS_2 = re.compile(
    r'^(\d+)(D)?([S-](\d+))?(/?DV(\d+)(D)?(-(\d))?)([\*\.](\d+))?$'
)


# --- HELPER CLASSES ---


class _DecomposedTCR(_DecomposedGene):
    @property
    def valid(self) -> bool:
        gene_str = self.compile(allele=False)

        # Is the gene valid?
        if not gene_str in self.allele_dict:
            return False
        
        # If the gene exists and an allele is specified, check if it exists
        if self.allele_num:
            return self.allele_num in self.allele_dict[gene_str]
        
        # If no allele specified, return true as the gene exists
        return True


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


class _DecomposedHomoSapiensTCR(_DecomposedTCR):
    allele_dict = _TCR_ALLELES_HOMOSAPIENS

    def __init__(
        self,
        base: Optional[str],
        num1: Optional[str],
        num2: Optional[str],
        p: bool,
        or92: bool,
        d_designation: Optional[str],
        allele_num: Optional[str]
    ) -> None:
        self.base = base
        self.num1 = num1
        self.num2 = num2
        self.p = p
        self.or92 = or92
        self.d_designation = d_designation
        self.allele_num = allele_num
    

    def compile(self, allele: bool = True) -> str:
        compiled = 'TR' + self.base + self.num1

        if self.num2:
            compiled = compiled + '-' + self.num2

        if self.p:
            compiled = compiled + 'P'

        if self.d_designation:
            compiled = compiled + '/DV' + self.d_designation
        
        if self.or92:
            compiled = compiled + '/OR9-2'
        
        if allele and self.allele_num:
            compiled = compiled + '*' + self.allele_num

        return compiled


    def resolve(self) -> bool:
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


class _DecomposedMusMusculusTCR(_DecomposedTCR):
    allele_dict = _TCR_ALLELES_MUSMUSCULUS

    def __init__(
        self,
        base: Optional[str],
        num1: Optional[str],
        num1_d: bool,
        num2: Optional[str],
        dv_1: Optional[str],
        dv_1_d: bool,
        dv_2: Optional[str],
        allele_num: Optional[str]
    ) -> None:
        self.base = base
        self.num1 = num1
        self.num1_d = num1_d
        self.num2 = num2
        self.dv_1 = dv_1
        self.dv_1_d = dv_1_d
        self.dv_2 = dv_2
        self.allele_num = allele_num
    

    def compile(self, allele: bool = True) -> str:
        compiled = 'TR' + self.base + self.num1

        if self.num1_d:
            compiled = compiled + 'D'

        if self.num2:
            compiled = compiled + '-' + self.num2

        if self.dv_1:
            compiled = compiled + '/DV' + self.dv_1

        if self.dv_1_d:
            compiled = compiled + 'D'

        if self.dv_2:
            compiled = compiled + '-' + self.dv_2
        
        if allele and self.allele_num:
            compiled = compiled + '*' + self.allele_num

        return compiled


    def resolve(self) -> None:
        # Remove any leading zeros in gene numbers
        self.num1 = None if self.num1 is None else self.num1.lstrip('0')
        self.num2 = None if self.num2 is None else self.num2.lstrip('0')

        # Resolve TRD designation if appropriate
        if self.dv_1 is None:
            try:
                d_designation = _TCR_D_DESIGNATIONS_MUSMUSCULUS[
                    self.compile(allele=False)
                ]
                self.dv_1 = d_designation['dv_1']
                self.dv_1_d = d_designation['dv_1_d']
                self.dv_2 = d_designation['dv_2']

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
    decomp_tcr = _DecomposedHomoSapiensTCR(
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


def _standardise_musmusculus(gene_name: str) -> str:
    # Take note of initial input for reference
    original_input = gene_name

    # Clean whitespace, remove known pollutors
    gene_name = ''.join(gene_name.split())
    gene_name = gene_name.replace('&nbsp;','')

    # Parse attempt 1
    if m := PARSE_RE_MUSMUSCULUS_1.match(gene_name): # ^TR([AB][VJ])(\d+)(D)?(-(\d+))?(/DV(\d+)(D)?(-(\d))?)?(\*(.+))?$
        base = m.group(1)
        num1 = m.group(2)
        num1_d = True if m.group(3) else False
        num2 = m.group(5)
        dv_1 = m.group(8)
        dv_1_d = True if m.group(9) else False
        dv_2 = m.group(11)
        allele_num = m.group(13)

    # Parse attempt 2
    elif m := PARSE_RE_MUSMUSCULUS_2.match(gene_name): # ^(\d+)(D)?([S-](\d+))?(/?DV(\d+)(D)?(-(\d))?)([\*\.](\d+))?$
        base = 'AV'
        num1 = m.group(1)
        num1_d = True if m.group(2) else False
        num2 = m.group(4)
        dv_1 = m.group(6)
        dv_1_d = True if m.group(7) else False
        dv_2 = m.group(9)
        allele_num = m.group(11)

    # Could not parse
    else:
        _warn_failure(original_input, gene_name, 'Mus musculus')
        return None

    # Build DecomposedTcr object
    decomp_tcr = _DecomposedMusMusculusTCR(
        base=base,
        num1=num1,
        num1_d=num1_d,
        num2=num2,
        dv_1=dv_1,
        dv_1_d=dv_1_d,
        dv_2=dv_2,
        allele_num=allele_num
    )

    # Try resolving, and return None on failure
    if not decomp_tcr.resolve():
        _warn_failure(original_input, decomp_tcr.compile(), 'Mus musculus')
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
    'HomoSapiens': _standardise_homosapiens,
    'MusMusculus': _standardise_musmusculus
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