'''
Gene standardiser classes
'''


from abc import ABC, abstractmethod
import re
from .._resources import *
from typing import List, Optional


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
    def valid(self) -> bool:
        '''
        Returns True if the gene's name is a real, valid and existing one.
        Else returns False.
        '''

    
    @abstractmethod
    def compile(self) -> str:
        '''
        Compile a complete string representation of the gene.
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
            r'^([A-Z0-9\-\.\(\)\/]+)(\*([0-9]+))?',
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
        if self.valid(enforce_functional=False):
            return # No resolution necessary

        # If a synonym, correct to currently approved name
        if self.syn_dict and self.gene in self.syn_dict:
            self.gene = self.syn_dict[self.gene]
            return # Done

        # Fix common errors
        self.gene = re.sub(r'(?<!TR)(?<!\/)DV', '/DV', self.gene)
        self.gene = re.sub(r'(?<!\d)0', '', self.gene)
        self.gene = self.gene.replace('TCR', 'TR')


    def valid(self, enforce_functional: bool) -> bool:
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
    

    def compile(self, precision: str) -> str:
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
            self.is_g = False
            self.is_p = False
            return

        parse_attempt = re.match(
            r'^([A-Z0-9\-\.\:\/]+)(\*([\d:]+G?P?)[LSCAQN]?)?',
            gene
        )
        
        if parse_attempt:
            self.gene = parse_attempt.group(1)
            self.allele_designation =\
                [] if parse_attempt.group(3) is None\
                else parse_attempt.group(3).split(':')
            
            if self.allele_designation:
                self.is_g = self.allele_designation[-1].endswith('G')
                self.is_p = self.allele_designation[-1].endswith('P')
            else:
                self.is_g = self.is_p = False
            
            return
        
        self.gene = gene
        self.allele_designation = None
        self.is_g = False
        self.is_p = False
    

    def resolve_errors(self) -> None:
        if self.valid():
            return # No resolution needed

        # If a synonym, correct to currently approved name
        if self.gene in HOMOSAPIENS_MHC_SYNONYMS:
            self.gene = HOMOSAPIENS_MHC_SYNONYMS[self.gene]
            if self.valid():
                return # Done
        
        # Add HLA prefix if necessary
        if not self.gene.startswith('HLA-'):
            self.gene = 'HLA-' + self.gene
            if self.valid():
                return # Done
        
        # Handle for forgotten asterisk
        if not self.allele_designation:
            m = re.match(r'^(HLA-[A-Z]+)(\d+)$', self.gene)
            if m:
                self.gene = m.group(1)
                self.allele_designation = [f'{int(m.group(2)):02}']


    def valid(self) -> bool:
        # Is the gene B2M?
        if self.gene == 'B2M' and not self.allele_designation:
            return True

        # Is the gene valid?
        if not self.gene in HOMOSAPIENS_MHC:
            return False

        # Verify allele designators up to the level of the protein (or G/P)
        allele_designation = self.allele_designation.copy()
        if not (self.is_g or self.is_p):
            allele_designation = allele_designation[:2]
        current_root = HOMOSAPIENS_MHC[self.gene]

        while len(allele_designation) > 0:
            try:
                current_root = current_root[allele_designation.pop(0)]
            except(KeyError):
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
    

    def compile(self, precision) -> str:
        if self.allele_designation:
            if precision == 'allele':
                return f'{self.gene}*{":".join(self.allele_designation)}'
            
            if precision == 'protein':
                prot_designation = self.allele_designation[:2]
                further_designation =\
                    None if len(self.allele_designation) <= 2\
                    else ':'+':'.join(self.allele_designation[2:])
                
                return (
                    f'{self.gene}*{":".join(prot_designation)}',
                    further_designation
                )

        return self.gene