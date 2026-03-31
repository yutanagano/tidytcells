from typing import Optional
from tidytcells._utils.alignment import get_compatible_symbols
from tidytcells._resources import SUPPORTED_RECEPTOR_SPECIES_AND_THEIR_AA_SEQUENCES


class ReceptorGene:
    '''
    A wrapper object for the receptor gene.

    If standardization was successful, this object provides access to the standardized allele/gene/subgroup and other properties.
    When failed, the error message(s) and attempted partially standardized gene symbol can be retrieved.
    '''
    def __init__(self, original_input, error, gene_name=None, allele_designation=None, subgroup_name=None, species=None):
        self._original_input = original_input
        self._error = error
        self._gene_name = gene_name
        self._allele_designation = allele_designation
        self._subgroup_name = subgroup_name
        self._species = species

        self._highest_precision_symbol = None

        if self._gene_name is not None and self._allele_designation is not None:
            self._highest_precision_symbol = f"{self._gene_name}*{self._allele_designation}"
        elif self._gene_name is not None:
            self._highest_precision_symbol = self._gene_name
        elif self._subgroup_name is not None:
            self._highest_precision_symbol = self._subgroup_name

    def __str__(self):
        str_repr = self.symbol

        if str_repr is not None:
            return str_repr
        else:
            return ""

    @property
    def original_input(self) -> Optional[str]:
        '''The original input symbol.'''
        return self._original_input

    @property
    def error(self) -> Optional[str]:
        '''The error message, only if standardization failed, otherwise ``None``.'''
        return self._error

    @property
    def is_standardized(self) -> bool:
        '''``True`` if the standardization was successful, ``False`` otherwise.'''
        return self.error is None

    @property
    def attempted_fix(self) -> Optional[str]:
        '''The best attempt at fixing the input symbol, only of standardization failed, if the standardization was a success this returns ``None``.'''
        if not self.is_standardized:
            return self._highest_precision_symbol

    @property
    def symbol(self) -> Optional[str]:
        '''The allele, gene or subgroup (whichever is most precise) if standardization was successful, otherwise ``None``.'''
        if self.is_standardized:
            return self._highest_precision_symbol

    @property
    def allele(self) -> Optional[str]:
        '''The allele name, if standardization was successful and allele-level information is available, otherwise ``None``.'''
        if self.is_standardized and self._allele_designation is not None and self._gene_name is not None:
            return f"{self._gene_name}*{self._allele_designation}"

    @property
    def gene(self) -> Optional[str]:
        '''The gene name, if standardization was successful and gene-level information is available, otherwise ``None``.'''
        if self.is_standardized:
            return self._gene_name

    @property
    def subgroup(self) -> Optional[str]:
        '''The subgroup name, if standardization was successful, otherwise ``None``.'''
        if self.is_standardized:
            return self._subgroup_name

    @property
    def locus(self) -> Optional[str]:
        '''
        The locus of the gene.
        This is typically the three-letter code ('TRA', 'TRB', 'TRG', 'TRD', 'IGH', 'IGL', 'IGK'),
        but for TRAV/DV genes, 'TRA/D' is returned.
        '''
        if self.is_standardized:
            locus = self.symbol[0:3]
            if "/D" in self.symbol:
                locus += "/D"
            return locus

    @property
    def receptor_type(self):
        ''''TR' for T cell receptor genes, or 'IG' for antibody genes if standardization was successful, otherwise ``None``.'''
        if self.is_standardized:
            return self.symbol[0:2]

    @property
    def gene_type(self) -> Optional[str]:
        '''The gene type ('V', 'D' or 'J'), if standardization was successful, otherwise ``None``.'''
        if self.is_standardized:
            return self.symbol[3]

    @property
    def species(self) -> str:
        '''The species used to validate the gene name.'''
        return self._species

    def get_all_alleles(self, enforce_functional=True):
        '''
        Get all alleles related to the standardized symbol

        :param enforce_functional:
            If ``True``, only functional alleles are returned
        :type enforce_functional:
            bool
        :return:
            A list of allele names
        :rtype:
            list
        '''
        if self.is_standardized:
            aa_dict = SUPPORTED_RECEPTOR_SPECIES_AND_THEIR_AA_SEQUENCES[self.receptor_type][self.species]

            return get_compatible_symbols(self.symbol, aa_dict, self.gene_type, self.locus, enforce_functional)

    def get_aa_sequences(self, sequence_type="ALL", enforce_functional=True):
        '''
        Get amino acid sequence information related to the alleles of the standardized symbol

        :param sequence_type:
            Which sequence to return. This can be:
            - For V genes: 'FR1', 'FR2', 'FR3', 'CDR1', 'CDR2', 'V-REGION'
            - For D genes: 'D-REGION'
            - For J genes: 'J-REGION', 'J-MOTIF'
            - Or 'ALL' to return all available sequences
        :type sequence_type:
            str
        :param enforce_functional:
            If ``True``, only information for functional alleles is returned
        :type enforce_functional:
            bool
        :return:
            A dictionary with allele names as keys and sequences as values
            When sequence_type is 'ALL', the result is a nested dictionary with allele names
            as outer keys, sequence types as inner keys, and sequences as inner values.
        :rtype:
            dict
        '''
        sequence_type = sequence_type.upper()
        sequence_type = sequence_type + "-IMGT" if sequence_type in {"FR1", "FR2", "FR3", "CDR1", "CDR2"} else sequence_type

        if self.is_standardized:
            aa_dict = SUPPORTED_RECEPTOR_SPECIES_AND_THEIR_AA_SEQUENCES[self.receptor_type][self.species]

            alleles_of_interest = self.get_all_alleles(enforce_functional)

            if sequence_type == "ALL":
                return {allele: aa_dict[allele] for allele in alleles_of_interest}
            else:
                return {allele: aa_dict[allele][sequence_type] if sequence_type in aa_dict[allele] else None
                        for allele in alleles_of_interest}
