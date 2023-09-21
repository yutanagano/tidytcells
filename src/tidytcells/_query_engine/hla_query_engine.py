from typing import FrozenSet
import warnings

from tidytcells._query_engine import QueryEngine
from tidytcells._resources import VALID_HOMOSAPIENS_MH


class HlaQueryEngine(QueryEngine):
    @classmethod
    def query(cls, precision: str, functionality: str) -> FrozenSet[str]:
        if precision == "allele":
            warnings.warn(
                "tidytcells is not fully aware of all HLA alleles, and the "
                "highest resolution it can provide is up to the level of the "
                "protein (two allele designations)."
            )

        query_results = []

        for gene_symbol, allele_dictionary in VALID_HOMOSAPIENS_MH.items():
            if precision == "gene":
                query_results.append(gene_symbol)
                continue

            for first_allele_designation in allele_dictionary:
                for second_allele_designation in allele_dictionary[
                    first_allele_designation
                ]:
                    if cls._is_g_or_p_group(second_allele_designation):
                        continue

                    query_results.append(
                        gene_symbol
                        + "*"
                        + first_allele_designation
                        + ":"
                        + second_allele_designation
                    )

        return frozenset(query_results)

    @classmethod
    def _is_g_or_p_group(cls, second_allele_designation: str) -> bool:
        return not second_allele_designation.isdigit()
