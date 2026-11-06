from abc import ABC
from abc import abstractmethod
from typing import Dict, Optional, Tuple, Set

from tidytcells.result import ReceptorGene



class ReceptorGeneSymbolStandardizer(ABC):

    def _is_synonym(self) -> bool:
        return self._gene_name in self._synonym_dictionary


    def __init__(self, symbol: str, enforce_functional: bool, allow_subgroup: bool,
                 species: str, germline_database: str,
                 valid_gene_dictionary: dict, synonym_dictionary: dict) -> None:
        self.original_symbol = symbol
        self.enforce_functional = enforce_functional

        self.species = species
        self.germline_database = germline_database


        self._valid_gene_dictionary = valid_gene_dictionary
        self._synonym_dictionary = {key: value for key, value in synonym_dictionary.items() if value in valid_gene_dictionary}
        self._valid_subgroups = {key.split("-")[0] for key in valid_gene_dictionary}


        self._gene_name, self._allele_designation = self._parse_symbol(symbol)
        self._allow_subgroup = self._allele_designation is None and allow_subgroup
        self._resolve_gene_name()
        self._compile_result()

    @abstractmethod
    def _parse_symbol(self, symbol: str) -> Tuple[Optional[str], Optional[str]]:
        pass

    @abstractmethod
    def _resolve_gene_name(self) -> None:
        pass

    def _has_valid_gene_name(self) -> bool:
        if self._gene_name in self._valid_gene_dictionary:
            return True

        if self._allow_subgroup and self._gene_name in self._valid_subgroups:
            return True

        return False

    def get_reason_why_invalid(self) -> Optional[str]:
        if not self._gene_name in self._valid_gene_dictionary:
            if self._gene_name in self._valid_subgroups:
                if self._allow_subgroup:
                    return None
                else:
                    return "Symbol is a subgroup (not a gene)"

            return "Unrecognized gene name"

        if self._allele_designation:
            allele_valid = (
                self._allele_designation in self._valid_gene_dictionary[self._gene_name]
            )

            if not allele_valid:
                return "Nonexistent allele for recognized gene"

            if (
                self.enforce_functional
                and self._valid_gene_dictionary[self._gene_name][self._allele_designation]
                != "F"
            ):
                return "Nonfunctional allele"

            return None

        if (
            self.enforce_functional
            and not "F" in self._valid_gene_dictionary[self._gene_name].values()
        ):
            return "Gene has no functional alleles"

        return None

    def _compile_result(self):
        gene_name = self._gene_name if self._gene_name in self._valid_gene_dictionary else None
        subgroup_name = self._gene_name.split("-")[0]
        subgroup_name = subgroup_name if subgroup_name in self._valid_subgroups else None

        self.result =  ReceptorGene(original_input=self.original_symbol,
                                    error=self.get_reason_why_invalid(),
                                    allele_designation=self._allele_designation,
                                    gene_name=gene_name,
                                    subgroup_name=subgroup_name,
                                    species=self.species,
                                    germline_database=self.germline_database)
