from typing import Optional


class MhGene:
    '''
    A wrapper object for the MHC gene.

    If standardization was successful, this object provides access to the standardized allele/gene and other properties.
    When failed, the error message(s) and attempted partially standardized gene symbol can be retrieved.
    '''
    def __init__(self, original_input, error, gene_name=None, allele_designation=None, species=None):
        self._original_input = original_input
        self._error = error
        self._gene_name = gene_name
        self._allele_designation = allele_designation if allele_designation is not None and len(allele_designation) > 0 else None
        self._species = species

        self._highest_precision_symbol = self._gene_name

        if self._gene_name is not None and self._allele_designation is not None:
            self._highest_precision_symbol = f'{self._gene_name}*{":".join(self._allele_designation)}'

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
        '''
        The best attempt at fixing the input symbol, only of standardization failed,
        if the standardization was a success this returns ``None``.
        '''
        if not self.is_standardized:
            return self._highest_precision_symbol

    @property
    def symbol(self) -> Optional[str]:
        '''The allele or gene (whichever is most precise) if standardization was successful, otherwise ``None``.'''
        if self.is_standardized:
            return self._highest_precision_symbol

    @property
    def allele(self) -> Optional[str]:
        '''The allele name, if standardization was successful and allele-level information is available, otherwise ``None``.'''
        if self.is_standardized and self._allele_designation is not None and self._gene_name is not None:
            return f'{self._gene_name}*{":".join(self._allele_designation)}'

    @property
    def gene(self) -> Optional[str]:
        '''The gene name, if standardization was successful, otherwise ``None``.'''
        if self.is_standardized and self._gene_name is not None:
            return self._gene_name

    @property
    def species(self) -> str:
        '''The species used to validate the gene name.'''
        return self._species


class HLAGene(MhGene):
    '''
    A wrapper object for the HLA gene.

    If standardization was successful, this object provides access to the standardized allele/protein/gene and other properties.
    When failed, the error message(s) and attempted partially standardized gene symbol can be retrieved.
    '''
    def __init__(self, original_input, error, gene_name=None, allele_designation=None):
        super().__init__(original_input, error, gene_name, allele_designation, species="homosapiens")

    @property
    def protein(self) -> Optional[str]:
        '''The protein name, if standardization was successful and protein-level information is available, otherwise ``None``.'''
        if self.is_standardized and self._allele_designation is not None and self._gene_name is not None:
            return f'{self._gene_name}*{":".join(self._allele_designation[:2])}'
