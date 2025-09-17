from abc import ABC, abstractmethod
from typing import Optional, Dict


class StandardizedJunction(ABC):
    """
    Abstract base standardizer class.
    """

    @abstractmethod
    def __init__(self, symbol: str) -> None:
        pass

    @property
    @abstractmethod
    def _synonym_dictionary(self) -> Dict[str, Dict]:
        pass

    @abstractmethod
    def get_reason_why_invalid(self) -> Optional[str]: # optional: enforce_complete / no reconstruction etc
        """
        If the CDR3 cannot be standardized (it is invalid), this method returns a string outlining the reason why (incomplete on the left side, right side, etc).
        Returns None if standardisation was successful.
        """

    @abstractmethod
    def compile(self, region: str = "junction") -> str:
        """
        Compile a complete string representation of the gene.
        The argument given to precision will determine the amount of specificity given in the compiled string.
        """
