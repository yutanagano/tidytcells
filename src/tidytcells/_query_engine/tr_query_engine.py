from abc import abstractmethod
from typing import Dict, FrozenSet

from tidytcells._query_engine import QueryEngine


class TrQueryEngine(QueryEngine):
    @property
    @abstractmethod
    def _valid_tr_dictionary(self) -> Dict[str, Dict[int, str]]:
        pass

    @classmethod
    def query(cls, precision: str, functionality: str) -> FrozenSet[str]:
        query_results = []

        for gene_symbol, allele_dictionary in cls._valid_tr_dictionary.items():
            if precision == "gene":
                if cls._gene_matches_functionality_requirements(
                    allele_dictionary, functionality
                ):
                    query_results.append(gene_symbol)
                continue

            for allele_designation, allele_functionality in allele_dictionary.items():
                if cls._allele_matches_functionality_requirements(
                    allele_functionality, functionality
                ):
                    query_results.append(gene_symbol + "*" + allele_designation)

        return frozenset(query_results)

    @classmethod
    def _gene_matches_functionality_requirements(
        cls, allele_dictionary: dict, functionality_setting: str
    ) -> bool:
        gene_alleles_matching_functionality_requirements = [
            cls._allele_matches_functionality_requirements(
                allele_functionality, functionality_setting
            )
            for allele_functionality in allele_dictionary.values()
        ]
        return any(gene_alleles_matching_functionality_requirements)

    @classmethod
    def _allele_matches_functionality_requirements(
        cls, allele_functionality: str, functionality_setting: str
    ) -> bool:
        if functionality_setting == "any":
            return True
        if functionality_setting == allele_functionality:
            return True
        if functionality_setting == "NF" and allele_functionality in ("P", "ORF"):
            return True
        return False
