from abc import abstractmethod
from typing import Dict, FrozenSet
from tidytcells._query_engine import QueryEngine
from tidytcells._resources import HOMOSAPIENS_MH_VALID_IPD, HOMOSAPIENS_MH_VALID_PROTEIN_MRO
import warnings


class HlaQueryEngine(QueryEngine):
    @property
    @abstractmethod
    def _valid_symbols(self) -> Dict[str, Dict]:
        pass

    @classmethod
    def query(cls, precision: str, functionality: str) -> FrozenSet[str]:
        if precision == "allele":
            warnings.warn(
                "tidytcells is not fully aware of all HLA alleles, and the "
                "highest resolution it can provide is up to the level of the "
                "protein (two allele designations)."
            )

        query_results = []

        for gene_symbol, allele_dictionary in cls._valid_symbols.items():
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


class ImgtHlaQueryEngine(HlaQueryEngine):
    _valid_symbols = HOMOSAPIENS_MH_VALID_IPD


class MroHlaQueryEngine(HlaQueryEngine):
    _valid_symbols = HOMOSAPIENS_MH_VALID_PROTEIN_MRO