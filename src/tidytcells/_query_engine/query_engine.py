from abc import ABC, abstractmethod
from typing import FrozenSet


class QueryEngine(ABC):
    @classmethod
    @abstractmethod
    def query(cls, precision: str, functionality: str) -> FrozenSet[str]:
        """
        Return the set of all unique genes/alleles of a particular type and species to the level of precision specified.
        The functionality argument specifies if query results should be filtered based on gene/allele functonality.
        """
