from typing import Optional

class Junction:
    '''
    A wrapper object for the CDR3 junction sequence.

    If standardization was successful, this object provides access to the standardized junction and CDR3.
    When failed, the error message(s) and attempted partially standardized CDR3 can be retrieved.
    '''

    def __init__(self, original_input, error, corrected_junction=None, species=None):
        self._original_input = original_input
        self._error = error
        self._corrected_junction = corrected_junction
        self._species = species

    def __str__(self):
        str_repr = self.junction

        if str_repr is not None:
            return str_repr
        else:
            return ""

    @property
    def original_input(self) -> Optional[str]:
        '''The original input sequence.'''
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
        The best attempt at fixing the input sequence, only of standardization failed,
        if the standardization was a success this returns ``None``.
        '''
        if not self.is_standardized:
            return self._corrected_junction

    @property
    def junction(self) -> Optional[str]:
        '''The IMGT-junction, including conserved leading C and trailing F / W / C if the standardization was successful, otherwise ``None``.'''
        if self.is_standardized:
            return self._corrected_junction

    @property
    def cdr3(self) -> Optional[str]:
        '''The IMGT-CDR3, excluding conserved leading C and trailing F / W / C if the standardization was successful, otherwise ``None``.'''
        if self.is_standardized:
            if self._corrected_junction is not None and len(self._corrected_junction) > 2:
                return self._corrected_junction[1:-1]

    @property
    def species(self) -> str:
        '''The species used for the gene lookup to validate the CDR3 junction.'''
        return self._species
