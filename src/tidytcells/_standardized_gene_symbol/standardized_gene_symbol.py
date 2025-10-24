from abc import ABC, abstractmethod
from typing import Optional


class StandardizedSymbol(ABC):
    """
    Abstract base standardizer class.
    """

    @abstractmethod
    def __init__(self, symbol: str, enforce_functional: bool, allow_subgroup: bool) -> None:
        pass

    @abstractmethod
    def get_reason_why_invalid(self, enforce_functional: bool = False) -> Optional[str]:
        """
        If the gene cannot be standardized (it is invalid), this method returns a string outlining the reason why (nonexistent gene, not functional, etc.).
        Returns None if standardisation was successful.
        If enforce_functional is set to True, then symbols of valid but non-functional genes will be rejected.
        """
