from abc import abstractmethod
import itertools
import re
from typing import Tuple, Optional
from typing import Dict, Optional, List, Set
from tidytcells import _utils
from tidytcells._standardized_gene_symbol import ReceptorGeneSymbolStandardizer


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


class TrSymbolStandardizer(ReceptorGeneSymbolStandardizer):

    def _parse_symbol(self, tr_symbol: str) -> Tuple[Optional[str], Optional[str]]:
        cleaned_tr_symbol = _utils.clean_and_uppercase(tr_symbol)
        parsed_tr_symbol = TrSymbolParser(cleaned_tr_symbol)

        return parsed_tr_symbol.gene_name, parsed_tr_symbol.allele_designation

    def _resolve_gene_name(self) -> None:
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

        self._try_resolving_trdv_designation_from_trav_info()
        if self._has_valid_gene_name():
            return

        self._try_resolving_trav_designation_from_trdv_info()
        if self._has_valid_gene_name():
            return

        self._try_removing_dash1()
        if self._has_valid_gene_name():
            return

    def _fix_common_errors_in_tr_gene_name(self) -> None:
        self._gene_name = self._gene_name.replace("TCR", "TR")
        self._gene_name = self._gene_name.replace("S", "-")
        self._gene_name = self._gene_name.replace(".", "-")
        self._gene_name = re.sub(r"(?<!TR)(?<!\/)-?DV", "/DV", self._gene_name)
        self._gene_name = re.sub(r"(?<!\/)-?OR", "/OR", self._gene_name)
        self._gene_name = re.sub(r"(?<!\d)0+", "", self._gene_name)

    def _try_resolving_trdv_designation_from_trav_info(self) -> None:
        if self._gene_name.startswith("TRAV") and "DV" not in self._gene_name:
            if "/" in self._gene_name:
                split_gene_name = self._gene_name.split("/")
                dv_segment = "DV" + split_gene_name.pop()
                self._gene_name = "/".join([*split_gene_name, dv_segment])
            else:
                for valid_gene_name in self._valid_gene_dictionary:
                    if valid_gene_name.startswith(self._gene_name + "/DV"):
                        self._gene_name = valid_gene_name
                        return

    def _try_resolving_trav_designation_from_trdv_info(self) -> None:
        if "DV" in self._gene_name:
            if self._gene_name.startswith("TRDV"):
                for valid_gene in self._valid_gene_dictionary:
                    dv_segment = self._gene_name[2:]
                    if re.match(rf"^TRAV\d+(-\d)?\/{dv_segment}$", valid_gene):
                        self._gene_name = valid_gene
            else:
                parse_attempt = re.match(r"^TR([\d-]+)\/(DV[\d-]+)$", self._gene_name)
                if parse_attempt:
                    self._gene_name = (
                        f"TRAV{parse_attempt.group(1)}/{parse_attempt.group(2)}"
                    )

    def _try_removing_dash1(self):
        orig_gene_name = self._gene_name

        if "-1" in self._gene_name:
            all_gene_nums = [
                (m.group(0), m.start(0), m.end(0))
                for m in re.finditer(r"\d+(-\d+)?", self._gene_name)
            ]

            dash1_candidates = []
            for numstr, start_idx, end_idx in all_gene_nums:
                fm = re.fullmatch(r"(\d+)(-1)?", numstr)
                if not fm:
                    continue
                dash1_candidates.append((fm.group(1), start_idx, end_idx))

            current_str_idx = 0
            without_dash1 = ""

            for numstr, start_idx, end_idx in dash1_candidates:
                without_dash1 += self._gene_name[current_str_idx:start_idx]
                without_dash1 += numstr
                current_str_idx = end_idx

            without_dash1 += self._gene_name[current_str_idx:]
            self._gene_name = without_dash1

            # only keep without_dash1 if this is a valid *gene* (don't convert from gene to subgroup)
            if self._gene_name in self._valid_gene_dictionary:
                return

            self._try_resolving_trdv_designation_from_trav_info()
            if self._gene_name in self._valid_gene_dictionary:
                return

        self._gene_name = orig_gene_name
