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


with resource_stream(__name__, 'resources/mhc_alleles_homosapiens.json') as s:
    _MHC_ALLELES_HOMOSAPIENS = json.load(s)
with resource_stream(__name__, 'resources/mhc_synonyms_homosapiens.json') as s:
    _MHC_SYNONYMS_HOMOSAPIENS = json.load(s)
PARSE_RE_HOMOSAPIENS = re.compile(
    r'^(HLA-)?([A-Za-z/]+\d??)(\*?([\d:]+(G)?(P)?)[LSCAQN]?)?$'
)

with resource_stream(__name__, 'resources/mhc_alleles_musmusculus.json') as s:
    _MHC_ALLELES_MUSMUSCULUS = json.load(s)
with resource_stream(__name__, 'resources/mhc_synonyms_musmusculus.json') as s:
    _MHC_SYNONYMS_MUSMUSCULUS = json.load(s)
PARSE_RE_MUSMUSCULUS = re.compile(r'^(M?H[12])?([A-Z0-9\.\/\']+)$')


# --- HELPER CLASSES ---


class _DecomposedHLA(_DecomposedGene):
    def __init__(
        self,
        gene: Optional[str],
        spec_fields: Optional[list],
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


    def resolve(self) -> bool:
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


class _DecomposedMusMusculusMHC(_DecomposedGene):
    def __init__(
        self,
        base: Optional[str],
        name: Optional[str]
    ) -> None:
        self.base = 'null' if base is None else base
        self.name = name


    @property
    def valid(self) -> bool:
        # Does the compiled gene name exist in the reference?
        if self.base in _MHC_ALLELES_MUSMUSCULUS:
            return self.name in _MHC_ALLELES_MUSMUSCULUS[self.base]
        
        return False
    

    def compile(self) -> Tuple[str, str]:
        prot_str = f'{self.base}-{self.name}'

        return (prot_str, None)


    def resolve(self) -> bool:
        if self.valid:
            return True
        
        # Fix deprecated names
        try:
            self.base, self.name =\
                _MHC_SYNONYMS_MUSMUSCULUS[self.base][self.name]

            # Deprecated name resolved
            return True
        except(KeyError):
            # Try sane modifications
            bases = ['H2', 'MH2', 'MH1']
            names = [self.name]

            if re.match(r'[A-Z]{2}', self.name):
                names.append(self.name[0])

            names += [name+'1' for name in names]

            for base, name in product(bases, names):
                self.base = base
                self.name = name

                if self.valid:
                    return True

        # Give up
        return False


# --- HELPER FUNCTIONS ---


def _standardise_homosapiens(gene_name: str) -> tuple[str]:
    species = 'Homo sapiens'

    # Take note of initial input for reference
    original_input = gene_name

    # Clean whitespace
    gene_name = ''.join(gene_name.split())

    # Parse attempt 1
    if gene_name == 'B2M':
        return ('B2M', None)

    # Parse attempt 2
    elif m := PARSE_RE_HOMOSAPIENS.match(gene_name): # ^(HLA-)?([A-Za-z/]+\d??)(\*?([\d:]+(G)?(P)?)[LSCAQN]?)?$
        # Extract gene name e.g. DRA1
        gene = m.group(2)
        # Extract the digits that specify which allele it is e.g. 01:01:01:01
        spec_fields = None if m.group(4) is None else m.group(4).split(':')
        g_group = True if m.group(5) else False
        p_group = True if m.group(6) else False

    # Could not parse
    else:
        _warn_failure(original_input, gene_name, species)
        return (None, None)

    # Build decomposed HLA object
    decomp_hla = _DecomposedHLA(
        gene=gene,
        spec_fields=spec_fields,
        g_group=g_group,
        p_group=p_group
    )

    # Try resolving, and return None on failure
    if not decomp_hla.resolve():
        _warn_failure(original_input, decomp_hla.compile(), species)
        return (None, None)
    
    return decomp_hla.compile()


def _standardise_musmusculus(gene_name: str) -> tuple[str]:
    species = 'Mus musculus'

    # Take note of initial input for reference
    original_input = gene_name

    # Clean whitespace
    gene_name = ''.join(gene_name.upper().split())
    gene_name = gene_name.replace('-', '')
    gene_name = gene_name.replace('ALPHA', 'A')
    gene_name = gene_name.replace('BETA', 'B')

    # Parse attempt 1
    if gene_name == 'B2M':
        return ('B2M', None)

    # Parse attempt 2
    elif m := PARSE_RE_MUSMUSCULUS.match(gene_name): # ^(M?H[12])?([A-Z0-9\.\/\']+)$
        base = m.group(1)
        gene = m.group(2)

    # Could not parse
    else:
        _warn_failure(original_input, gene_name, species)
        return (None, None)

    # Build decomposed HLA object
    decomp_mhc = _DecomposedMusMusculusMHC(
        base=base,
        name=gene
    )

    # Try resolving, and return None on failure
    if not decomp_mhc.resolve():
        _warn_failure(original_input, decomp_mhc.compile(), species)
        return (None, None)
    
    return decomp_mhc.compile()


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
    'HomoSapiens': _standardise_homosapiens,
    'MusMusculus': _standardise_musmusculus
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
        _warn_failure(gene_name, gene_name, species)
        return (None, None)

    # If the specified species is supported, attempt standardisation
    if species in SUPPORTED_SPECIES:
        return SUPPORTED_SPECIES[species](gene_name=gene_name)
    
    # Otherwise, don't touch it
    warn(
        f'Unsupported species: "{species}". '
        'Skipping MHC gene standardisation procedure...'
    )
    return (gene_name, None)


def get_chain(gene_name: str) -> str:
    '''
    Given a standardised MHC gene name, detect whether it codes for an alpha
    or a beta chain molecule.
    
    :param gene_name: Standardised MHC gene name (non-standardised values will
        be unrecognised)
    :type gene_name: str
    :return: ``'alpha'`` or ``'beta'`` if ``gene_name`` is recognised, else
        ``None``.
    :rtype: str or None

    '''

    if type(gene_name) == str:
        if re.match('HLA-([ABCEFG]|D[PQR]A)', gene_name):
            return 'alpha'
        
        if re.match('HLA-D[PQR]B|B2M', gene_name):
            return 'beta'

    warn(
        f'MHC gene unrecognised: "{gene_name}" '
        f'(type: {type(gene_name)}), cannot determine if alpha or beta.'
    )


def classify(gene_name: str) -> int:
    '''
    Given a standardised MHC gene name, detect whether it comprises a class I
    or II MHC receptor complex.
    
    :param gene_name: Standardised MHC gene name (non-standardised values will
        be unrecognised)
    :type gene_name: str
    :return: ``1`` or ``2`` if ``gene_name`` is recognised, else ``None``.
    :rtype: int or None

    '''

    if type(gene_name) == str:
        if re.match('HLA-[ABCEFG]|B2M', gene_name):
            return 1
        
        if re.match('HLA-D[PQR][AB]', gene_name):
            return 2

    warn(
        f'MHC gene unrecognised: "{gene_name}" '
        f'(type: {type(gene_name)}), cannot determine if class 1 or 2.'
    )