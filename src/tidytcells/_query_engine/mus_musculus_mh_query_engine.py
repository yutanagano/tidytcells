from typing import FrozenSet
import warnings

from tidytcells._query_engine import QueryEngine
from tidytcells._resources import MUSMUSCULUS_MHC


class MusMusculusMhQueryEngine(QueryEngine):
    @classmethod
    def query(cls, precision: str, functionality: str) -> FrozenSet[str]:
        if precision == "allele":
            warnings.warn(
                "tidytcells is not aware of Mus musculus MHC alleles at all, "
                "and can only provide up to the level of the gene."
            )

        return frozenset(MUSMUSCULUS_MHC)
