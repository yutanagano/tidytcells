from abc import ABC, abstractmethod
from typing import Optional


class StandardizedGeneSymbol(ABC):
    """
    Abstract base standardizer class.
    """

    @abstractmethod
    def __init__(self, gene_symbol: str) -> None:
        pass

    @abstractmethod
    def get_reason_why_invalid(self, enforce_functional: bool = False) -> Optional[str]:
        """
        If the gene cannot be standardized (it is invalid), this method returns a string outlining the reason why (nonexistent gene, not functional, etc.).
        Returns None if standardisation was successful.
        If enforce_functional is set to True, then symbols of valid but non-functional genes will be rejected.
        """

    @abstractmethod
    def compile(self, precision: str = "allele") -> str:
        """
        Compile a complete string representation of the gene.
        The argument given to precision will determine the amount of specificity given in the compiled string.
        """
