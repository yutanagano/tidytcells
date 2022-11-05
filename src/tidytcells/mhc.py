'Utility functions related to MHCs and MHC genes.'


import json
from pkg_resources import resource_stream
import re
from typing import Tuple, Union
from warnings import warn


# --- STATIC RESOURCES ---


with resource_stream(__name__, 'resources/mhc_alleles_homosapiens.json') as s:
    _MHC_ALLELES_HOMOSAPIENS = json.load(s)
with resource_stream(__name__, 'resources/mhc_synonyms_homosapiens.json') as s:
    _MHC_SYNONYMS_HOMOSAPIENS = json.load(s)
PARSE_RE_HOMOSAPIENS = re.compile(
    r'^(HLA-)?([A-Za-z/]+\d??)(\*?([\d:]+(G)?(P)?)[LSCAQN]?)?$'
)


# --- HELPER CLASSES ---


class DecomposedHla:
    def __init__(
        self,
        gene: Union[str, None],
        spec_fields: Union[list, None],
        g_group: bool,
        p_group: bool
    ) -> None:
        self.gene = gene
        self.spec_fields = [] if spec_fields is None else spec_fields
        self.g_group = g_group
        self.p_group = p_group


    @property
    def valid(self) -> bool:
        # Is the gene valid?
        if not self.gene in _MHC_ALLELES_HOMOSAPIENS:
            return False

        # If the gene exists and there are further specifier fields, walk down
        # allele tree up to the level of the protein (or the level of the G or
        # P group when appropriate) to see if the specifier field values are
        # valid
        spec_fields = self.spec_fields.copy()
        if not (self.g_group or self.p_group):
            spec_fields = spec_fields[:2]
        current_root = _MHC_ALLELES_HOMOSAPIENS[self.gene]

        while len(spec_fields) > 0:
            try:
                current_root = current_root[spec_fields.pop(0)]
            except(KeyError):
                # We exit the tree at some point, so specifier fields are
                # invalid
                return False

        # If there are specifier fields past the protein level, just make sure
        # they look like legitimate specifier field values
        if not (self.g_group or self.p_group) and len(self.spec_fields) > 2:
            further_specifiers = self.spec_fields[2:]

            if len(further_specifiers) > 2:
                return False

            for field in further_specifiers:
                if not field.isdigit():
                    return False
                
                if len(field) < 2:
                    return False
        
        # Valid gene with valid specifier fields, so valid!
        return True
    

    def compile(self) -> Tuple[str, str]:
        prot_str = 'HLA-' + self.gene
        
        if self.g_group or self.p_group:
            prot_str = prot_str + '*' + ':'.join(self.spec_fields)
            return (prot_str, None)

        if len(self.spec_fields) > 0:
            prot_str = prot_str + '*' + ':'.join(self.spec_fields[:2])
        
        if len(self.spec_fields) > 2:
            mut_str = ':' + ':'.join(self.spec_fields[2:])
        else:
            mut_str = None

        return (prot_str, mut_str)


    def resolve(self) -> None:
        if self.gene == 'Cw':
            self.gene = 'C'

        # If gene name looks deprecated, replace with modern counterpart
        if self.gene in _MHC_SYNONYMS_HOMOSAPIENS:
            self.gene = _MHC_SYNONYMS_HOMOSAPIENS[self.gene]

        # Clean spec fields
        self._resolve_spec_fields()

        # If now a valid HLA, return True
        if self.valid:
            return True

        # If still invalid, try adding or removing a 1 at the end of the gene
        # name to see if this makes it valid
        if self.gene.endswith('1'):
            self.gene = self.gene.rstrip('1')

            if self.valid:
                return True

            self.gene = self.gene + '1'
        else:
            self.gene = self.gene + '1'

            if self.valid:
                return True
            
            self.gene = self.gene.rstrip('1')

        return False


    def _resolve_spec_fields(self) -> None:
        # Clean up spec fields if any
        if len(self.spec_fields) > 0:
            # Colon is most likely missing between first two fields e.g.
            # HLA-A*0101 <-> HLA-A*01:01 so we should separate the four-digit
            # field in half
            if len(self.spec_fields[0]) == 4:
                field_0_original = self.spec_fields.pop(0)
                field_0 = field_0_original[:2]
                field_1 = field_0_original[2:]
                self.spec_fields = [field_0, field_1] + self.spec_fields
            
            # Add leading zeros to elements if necessary
            self.spec_fields = [
                '0'+f if len(f) == 1 else f for f in self.spec_fields
            ]


# --- HELPER FUNCTIONS ---


def _standardise_homosapiens(gene_str: str) -> str:
    # Take note of initial input for reference
    original_input = gene_str

    # Clean whitespace
    gene_str = ''.join(gene_str.split())

    # Parse attempt 1
    if gene_str == 'B2M':
        return ('B2M', None)

    # Parse attempt 1
    elif m := PARSE_RE_HOMOSAPIENS.match(gene_str):
        # Extract gene name e.g. DRA1
        gene = m.group(2)
        # Extract the digits that specify which allele it is e.g. 01:01:01:01
        spec_fields = None if m.group(4) is None else m.group(4).split(':')
        g_group = True if m.group(5) else False
        p_group = True if m.group(6) else False

    # Could not parse
    else:
        _warn_failure(original_input, gene_str, 'Homo sapiens')
        return (None, None)

    # Build decomposed HLA object
    decomp_hla = DecomposedHla(
        gene=gene,
        spec_fields=spec_fields,
        g_group=g_group,
        p_group=p_group
    )

    # Try resolving, and return None on failure
    if not decomp_hla.resolve():
        _warn_failure(original_input, decomp_hla.compile(), 'Homo sapiens')
        return (None, None)
    
    return decomp_hla.compile()


def _warn_failure(
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
    'HomoSapiens': _standardise_homosapiens
}


def standardise(gene_str: str, species: str) -> tuple:
    # If gene_str is not a string, skip and return None.
    if type(gene_str) != str:
        _warn_failure(gene_str, gene_str, species)
        return (None, None)

    # If the specified species is supported, attempt standardisation
    if species in SUPPORTED_SPECIES:
        return SUPPORTED_SPECIES[species](gene_str=gene_str)
    
    # Otherwise, don't touch it
    warn(
        f'Unsupported species: "{species}". '
        'Skipping MHC gene standardisation procedure...'
    )
    return (gene_str, None)


def get_chain(mhc_gene_name: str) -> str:
    '''
    Given an MHC gene name, classify it as MHC A or MHC B. NOTE: Currently only
    considers human MHCs, and ignores mouse.
    '''

    if type(mhc_gene_name) == str:
        if re.match('HLA-([ABCEFG]|D[PQR]A)', mhc_gene_name):
            return 'alpha'
        
        if re.match('HLA-D[PQR]B|B2M', mhc_gene_name):
            return 'beta'

    warn(
        f'MHC gene unrecognised: "{mhc_gene_name}" '
        f'(type: {type(mhc_gene_name)}), cannot determine if alpha or beta.'
    )


def classify(mhc_gene_name: str) -> int:
    '''
    Given an MHC gene name, classify it as MHC class 1 or 2. NOTE: Currently
    only considers HLA.
    '''

    if type(mhc_gene_name) == str:
        if re.match('HLA-[ABCEFG]|B2M', mhc_gene_name):
            return 1
        
        if re.match('HLA-D[PQR][AB]', mhc_gene_name):
            return 2

    warn(
        f'MHC gene unrecognised: "{mhc_gene_name}" '
        f'(type: {type(mhc_gene_name)}), cannot determine if class 1 or 2.'
    )