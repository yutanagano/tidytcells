'''
Gene standardiser classes
'''


from abc import ABC, abstractmethod
from itertools import product
import re
from .._resources import *


# --- STATIC RESOURCES ---


SUB_DV_RE = re.compile(r'(?<!TR)(?<!\/)DV')
SUB_OR_RE = re.compile(r'(?<!\/)OR')
SUB_ZERO_RE = re.compile(r'(?<!\d)0')


# --- STANDARDISER CLASSES ---


class GeneStandardiser(ABC):
    '''
    Abstract base standardiser class.
    '''


    @abstractmethod
    def __init__(self, gene: str) -> None:
        '''
        Init method for the standardiser objects. Should expect a whitespace-
        cleaned, uppercased string to be supplied to the 'gene' parameter.
        '''


    @abstractmethod
    def valid(self, enforce_functional: bool = False) -> bool:
        '''
        Returns True if the gene's name is a real, valid and existing one.
        Else returns False. If enforce_functional, additionally verify that
        the gene in question is functional.
        '''

    
    @abstractmethod
    def compile(self, precision: str = 'allele') -> str:
        '''
        Compile a complete string representation of the gene. The argument
        given to precision will determine the amount of specificity given
        in the compiled string.
        '''


class TCRStandardiser(GeneStandardiser):
    '''
    TCR standardiser base class.
    '''


    ref_dict = None
    syn_dict = None


    def __init__(self, gene: str) -> None:
        self.parse_gene_str(gene)
        self.resolve_errors()


    def parse_gene_str(self, gene: str) -> None:
        parse_attempt = re.match(
            r'^([A-Z0-9\-\.\(\)\/]+)(\*(\d+))?',
            gene
        )
        
        if parse_attempt:
            self.gene = parse_attempt.group(1)
            self.allele_designation =\
                None if parse_attempt.group(3) is None\
                else f'{int(parse_attempt.group(3)):02}'
            return
        
        self.gene = gene
        self.allele_designation = None


    def resolve_errors(self) -> None:
        if self.valid():
            return # No resolution necessary

        # If a synonym, correct to currently approved name
        if self.syn_dict and self.gene in self.syn_dict:
            self.gene = self.syn_dict[self.gene]
            return

        # Fix common errors
        self.gene = re.sub(r'(?<!TR)(?<!\/)DV', '/DV', self.gene)
        self.gene = re.sub(r'(?<!\d)0', '', self.gene)
        self.gene = self.gene.replace('TCR', 'TR')
        if self.valid():
            return
        
        # Make sure gene starts with 'TR'
        if not self.gene.startswith('TR'):
            self.gene = 'TR' + self.gene
            if self.valid():
                return
        
        # Resolve DV designation from AV if necessary
        if self.gene.startswith('TRAV') and 'DV' not in self.gene:
            for valid_gene in self.ref_dict:
                if valid_gene.startswith(self.gene + '/DV'):
                    self.gene = valid_gene
                    return
        
        # Resolve AV designation from DV is necessary
        if self.gene.startswith('TRDV'):
            for valid_gene in self.ref_dict:
                if re.match(
                    rf'^TRAV\d+(-\d)?\/{self.gene[2:]}$',
                    valid_gene
                ):
                    self.gene = valid_gene
                    return


    def valid(self, enforce_functional: bool = False) -> bool:
        # Is the gene valid?
        if not self.gene in self.ref_dict:
            return False
        
        # If the gene exists and an allele is specified, check if it exists
        if self.allele_designation:
            allele_valid = self.allele_designation in self.ref_dict[self.gene]

            if not enforce_functional:
                return allele_valid

            return self.ref_dict[self.gene][self.allele_designation] == 'F'
        
        # If enforce_functional, ensure there is at least one functional allele
        if enforce_functional:
            return 'F' in self.ref_dict[self.gene].values()
        
        # Otherwise gene is valid so return true
        return True
    

    def compile(self, precision: str = 'allele') -> str:
        if precision == 'allele' and self.allele_designation:
            return f'{self.gene}*{self.allele_designation}'

        return self.gene


class HomoSapiensTCRStandardiser(TCRStandardiser):
    ref_dict = HOMOSAPIENS_TCR
    syn_dict = HOMOSAPIENS_TCR_SYNONYMS


    def resolve_errors(self) -> None:
        super().resolve_errors()

        # Fix common errors
        self.gene = re.sub(r'(?<!\/)OR', '/OR', self.gene)


class MusMusculusTCRStandardiser(TCRStandardiser):
    ref_dict = MUSMUSCULUS_TCR


class HLAStandardiser(GeneStandardiser):
    def __init__(self, gene: str) -> None:
        self.parse_gene(gene)
        self.resolve_errors()

    
    def parse_gene(self, gene: str) -> None:
        if gene == 'B2M':
            self.gene = 'B2M'
            self.allele_designation = []
            return
        
        # If period between digits, replace with colon
        gene = re.sub(r'(?<=\d)\.(?=\d)', ':', gene)

        parse_attempt_1 = re.match(
            r'^((HLA-)?(D[PQ][AB]|DRB|TAP)\d)(\*?([\d:]+G?P?)[LSCAQN]?)?',
            gene
        )
        if parse_attempt_1:
            self.gene = parse_attempt_1.group(1)
            self.allele_designation =\
                [] if parse_attempt_1.group(5) is None\
                else parse_attempt_1.group(5).split(':')
            return

        parse_attempt_2 = re.match(
            r'^([A-Z0-9\-\.\:\/]+)(\*([\d:]+G?P?)[LSCAQN]?)?',
            gene
        )
        if parse_attempt_2:
            self.gene = parse_attempt_2.group(1)
            self.allele_designation =\
                [] if parse_attempt_2.group(3) is None\
                else parse_attempt_2.group(3).split(':')
            return
        
        self.gene = gene
        self.allele_designation = []
    

    @property
    def is_group(self) -> bool:
        if not self.allele_designation:
            return False
        
        return self.allele_designation[-1].endswith('G') or \
            self.allele_designation[-1].endswith('P')
    

    def resolve_errors(self) -> None:
        if self.valid():
            return # No resolution needed

        # If a synonym, correct to currently approved name
        if self.gene in HOMOSAPIENS_MHC_SYNONYMS:
            self.gene = HOMOSAPIENS_MHC_SYNONYMS[self.gene]
            if self.valid():
                return
        
        # Handle common errors
        if not self.gene.startswith('HLA-'):
            self.gene = 'HLA-' + self.gene
        self.gene = self.gene.replace('CW', 'C')
        if self.valid():
            return
        
        # Handle for forgotten asterisk
        if not self.allele_designation:
            m = re.match(r'^(HLA-[A-Z]+)([\d:]+G?P?)$', self.gene)
            if m:
                self.gene = m.group(1)
                self.allele_designation = m.group(2).split(':')
            if self.valid():
                return
        
        # Handle forgotten colon between first/second allele designator
        if self.allele_designation and len(self.allele_designation[0]) == 4:
            self.allele_designation = [
                self.allele_designation[0][:2],
                self.allele_designation[0][2:]
            ] + self.allele_designation[1:]
            if self.valid():
                return

        # Try different amounts of leading zeros in first 2 allele designators
        ads = [int(ad) for ad in self.allele_designation[:2]]
        ads_reformatted = [[f'{ad:02}', f'{ad:03}'] for ad in ads]
        for new_ads in product(*ads_reformatted):
            self.allele_designation =\
                list(new_ads) + self.allele_designation[2:]
            if self.valid():
                return


    def valid(self, enforce_functional: bool = False) -> bool:
        # Is the gene B2M?
        if self.gene == 'B2M' and not self.allele_designation:
            return True

        # Is the gene valid?
        if not self.gene in HOMOSAPIENS_MHC:
            return False

        # Verify allele designators up to the level of the protein (or G/P)
        allele_designation = self.allele_designation.copy()
        if not self.is_group:
            allele_designation = allele_designation[:2]
        current_root = HOMOSAPIENS_MHC[self.gene]

        while len(allele_designation) > 0:
            try:
                current_root = current_root[allele_designation.pop(0)]
            except(KeyError):
                return False

        # If there are designator fields past the protein level, just make sure
        # they look like legitimate designator field values
        if not self.is_group and len(self.allele_designation) > 2:
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
    

    def compile(self, precision) -> str:
        if self.allele_designation:
            if precision == 'allele':
                return f'{self.gene}*{":".join(self.allele_designation)}'
            
            if precision == 'protein':
                return f'{self.gene}*{":".join(self.allele_designation[:2])}'

        return self.gene


class MusMusculusMHCStandardiser(GeneStandardiser):
    def __init__(self, gene: str) -> None:
        self.parse_gene(gene)
        self.resolve_errors()
    

    def parse_gene(self, gene: str) -> None:
        parse_attempt = re.match(
            r'^([A-Z0-9\-\.\(\)\/]+)(\*(\d+))?',
            gene
        )

        if parse_attempt:
            self.gene = parse_attempt.group(1)
            self.allele_designation =\
                None if parse_attempt.group(3) is None\
                else f'{int(parse_attempt.group(3)):02}'
            return
        
        self.gene = gene
        self.allele_designation = None
    

    def resolve_errors(self) -> None:
        if self.valid():
            return # No resolution needed

        # If a synonym, correct to currently approved name
        if self.gene.replace('-', '') in MUSMUSCULUS_MHC_SYNONYMS:
            self.gene = MUSMUSCULUS_MHC_SYNONYMS[self.gene.replace('-', '')]
            if self.valid():
                return
    

    def valid(self, enforce_functional: bool = False) -> bool:
        return self.gene in MUSMUSCULUS_MHC


    def compile(self, precision: str = 'allele') -> str:
        if precision == 'allele' and self.allele_designation:
            return f'{self.gene}*{self.allele_designation}'

        return self.gene