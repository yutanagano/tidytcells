from abc import abstractmethod
import re
from typing import Dict, Optional

from tidytcells import _utils
from tidytcells._standardized_gene_symbol import StandardizedGeneSymbol


class TrSymbolParser:
    gene_name: str
    allele_designation: int

    def __init__(self, tr_symbol: str) -> None:
        parse_attempt = re.match(r"^([A-Z0-9\-\.\(\)\/]+)(\*(\d+))?", tr_symbol)

        if parse_attempt:
            self.gene_name = parse_attempt.group(1)
            self.allele_designation = (
                None
                if parse_attempt.group(3) is None
                else f"{int(parse_attempt.group(3)):02}"
            )
        else:
            self.gene_name = tr_symbol
            self.allele_designation = None


class StandardizedTrSymbol(StandardizedGeneSymbol):
    @property
    @abstractmethod
    def _synonym_dictionary(self) -> Dict[str, str]:
        pass

    @property
    @abstractmethod
    def _valid_tr_dictionary(self) -> Dict[str, Dict[int, str]]:
        pass

    def __init__(self, gene_symbol: str) -> None:
        self._parse_tr_symbol(gene_symbol)
        self._resolve_gene_name()

    def _parse_tr_symbol(self, tr_symbol: str) -> None:
        cleaned_tr_symbol = _utils.clean_and_uppercase(tr_symbol)
        parsed_tr_symbol = TrSymbolParser(cleaned_tr_symbol)
        self._gene_name = parsed_tr_symbol.gene_name
        self._allele_designation = parsed_tr_symbol.allele_designation

    def _resolve_gene_name(self, skip_dash1_section: bool = False) -> None:
        if self._has_valid_gene_name():
            return

        if self._is_synonym():
            self._gene_name = self._synonym_dictionary[self._gene_name]
            return

        self._fix_common_errors_in_tr_gene_name()
        if self._has_valid_gene_name():
            return

        if not self._gene_name.startswith("TR"):
            self._gene_name = "TR" + self._gene_name
            if self._has_valid_gene_name():
                return

        if self._gene_name.startswith("TRAV") and "DV" not in self._gene_name:
            self._try_resolving_trdv_designation_from_trav_info()
            if self._has_valid_gene_name():
                return

        if "DV" in self._gene_name:
            self._try_resolving_trav_designation_from_trdv_info()
            if self._has_valid_gene_name():
                return

        if not skip_dash1_section:
            original = self._gene_name
            self._try_adding_or_removing_dash1_to_gene_name()
            self._resolve_gene_name(skip_dash1_section=True)
            if self._has_valid_gene_name():
                return
            self._gene_name = original

    def _has_valid_gene_name(self) -> bool:
        return self._gene_name in self._valid_tr_dictionary

    def _is_synonym(self) -> bool:
        return self._gene_name in self._synonym_dictionary

    def _fix_common_errors_in_tr_gene_name(self) -> None:
        self._gene_name = self._gene_name.replace("TCR", "TR")
        self._gene_name = self._gene_name.replace("S", "-")
        self._gene_name = self._gene_name.replace(".", "-")
        self._gene_name = re.sub(r"(?<!TR)(?<!\/)DV", "/DV", self._gene_name)
        self._gene_name = re.sub(r"(?<!\/)OR", "/OR", self._gene_name)
        self._gene_name = re.sub(r"(?<!\d)0", "", self._gene_name)

    def _try_resolving_trdv_designation_from_trav_info(self) -> None:
        if "/" in self._gene_name:
            split_gene_name = self._gene_name.split("/")
            dv_segment = "DV" + split_gene_name.pop()
            self._gene_name = "/".join([*split_gene_name, dv_segment])
        else:
            for valid_gene_name in self._valid_tr_dictionary:
                if valid_gene_name.startswith(self._gene_name + "/DV"):
                    self._gene_name = valid_gene_name

    def _try_resolving_trav_designation_from_trdv_info(self) -> None:
        if self._gene_name.startswith("TRDV"):
            for valid_gene in self._valid_tr_dictionary:
                dv_segment = self._gene_name[2:]
                if re.match(rf"^TRAV\d+(-\d)?\/{dv_segment}$", valid_gene):
                    self._gene_name = valid_gene
        else:
            parse_attempt = re.match(r"^TR([\d-]+)\/(DV[\d-]+)$", self._gene_name)
            if parse_attempt:
                self._gene_name = (
                    f"TRAV{parse_attempt.group(1)}/{parse_attempt.group(2)}"
                )

    def _try_adding_or_removing_dash1_to_gene_name(self) -> None:
        if "-1" in self._gene_name:
            self._gene_name = self._gene_name.replace("-1", "")
        else:
            self._gene_name += "-1"

    def get_reason_why_invalid(self, enforce_functional: bool = False) -> Optional[str]:
        if not self._gene_name in self._valid_tr_dictionary:
            return "unrecognised gene name"

        if self._allele_designation:
            allele_valid = (
                self._allele_designation in self._valid_tr_dictionary[self._gene_name]
            )

            if not allele_valid:
                return "nonexistent allele for recognised gene"

            if (
                enforce_functional
                and self._valid_tr_dictionary[self._gene_name][self._allele_designation]
                != "F"
            ):
                return "nonfunctional allele"

            return None

        if (
            enforce_functional
            and not "F" in self._valid_tr_dictionary[self._gene_name].values()
        ):
            return "gene has no functional alleles"

        return None

    def compile(self, precision: str = "allele") -> str:
        if precision == "allele" and self._allele_designation is not None:
            return f"{self._gene_name}*{self._allele_designation}"

        return self._gene_name
