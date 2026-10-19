import re
import itertools
from typing import Tuple, Optional
from tidytcells import _utils
from tidytcells._standardized_gene_symbol import ReceptorGeneSymbolStandardizer


class IgSymbolParser:
    gene_name: str
    allele_designation: Optional[str]

    def __init__(self, ig_symbol: str, species: str, germline_database: str) -> None:
        if species == "homosapiens":
            gene_pattern = r"^([A-Z0-9\-\.\(\)\/ab]+)"
        elif species == "musmusculus":
            gene_pattern = r"^([A-Z0-9\-\.\(\)\/]+[abcdef]?)"
        else:
            gene_pattern = r"^([A-Z0-9\-\.\(\)\/]+)"

        if germline_database == "OGRDB":
            allele_pattern =  r"(\*(i?\d+))?"
        else:
            allele_pattern = r"(\*(\d+))?"

        pattern = gene_pattern + allele_pattern

        parse_attempt = re.match(pattern, ig_symbol)

        if parse_attempt:
            self.gene_name = parse_attempt.group(1)

            allele = parse_attempt.group(3)

            if allele is None:
                self.allele_designation = None
            elif allele.startswith("i"):
                self.allele_designation = f"i{int(allele[1:]):02}"
            else:
                self.allele_designation = f"{int(allele):02}"
        else:
            self.gene_name = ig_symbol
            self.allele_designation = None


class IgSymbolStandardizer(ReceptorGeneSymbolStandardizer):

    def _parse_symbol(self, ig_symbol: str) -> Tuple[Optional[str], Optional[str]]:
        cleaned_ig_symbol = self._safe_clean_ig_symbol(ig_symbol)
        parsed_ig_symbol = IgSymbolParser(cleaned_ig_symbol, species=self._species, germline_database=self._germline_database)

        return parsed_ig_symbol.gene_name, parsed_ig_symbol.allele_designation

    def _safe_clean_ig_symbol(self, ig_symbol: str) -> str:
        gene_part, star, allele_part = ig_symbol.partition("*")
        gene_part = _utils.clean_and_uppercase(gene_part)
        allele_part = _utils.clean_and_lowercase(allele_part)

        cleaned_ig_symbol = gene_part + star + allele_part

        # Deal with lowercase a/b in genes like IGHD1/OR15-1a*01
        if self._species == "homosapiens":
            match = re.search(r"(.*?OR15-\d)([AB])", cleaned_ig_symbol)
            if match:
                cleaned_ig_symbol = (
                    match.group(1)
                    + match.group(2).lower()
                    + cleaned_ig_symbol[match.end() :]
                )
        elif self._species == "musmusculus":
            match = re.search(r"(IGKV\d+-)([ABCDEF])", cleaned_ig_symbol)
            if match:
                cleaned_ig_symbol = (
                        match.group(1)
                        + match.group(2).lower()
                        + cleaned_ig_symbol[match.end():]
                )

        return cleaned_ig_symbol

    def _resolve_gene_name(self) -> None:
        if self._has_valid_gene_name():
            return

        if self._is_synonym():
            self._gene_name = self._synonym_dictionary[self._gene_name]
            return

        self._fix_common_errors_in_ig_gene_name()
        if self._has_valid_gene_name():
            return

        if not self._gene_name.startswith("IG"):
            self._gene_name = "IG" + self._gene_name
            if self._has_valid_gene_name():
                return

        self._try_removing_dash1()
        if self._has_valid_gene_name():
            return

    def _fix_common_errors_in_ig_gene_name(self) -> None:
        self._gene_name = self._gene_name.replace(".", "-")
        self._gene_name = re.sub(r"(?<!\/)-?OR", "/OR", self._gene_name)
        self._gene_name = re.sub(r"(?<!\d)0+", "", self._gene_name)

    def _try_removing_dash1(self):
        if "-1" not in self._gene_name:
            return

        original_gene_name = self._gene_name

        all_gene_nums = [
            (m.group(0), m.start(0), m.end(0))
            for m in re.finditer(r"\d+(-\d+)?", self._gene_name)
        ]

        dash1_candidates = []
        for numstr, start_idx, end_idx in all_gene_nums:
            fm = re.fullmatch(r"(\d+)(-1)", numstr)
            if not fm:
                continue
            dash1_candidates.append((fm.group(1), start_idx, end_idx))

        dash1_variants = []
        for comb in itertools.product(("keep", "remove"), repeat=len(dash1_candidates)):
            num_comb_zip = zip(dash1_candidates, comb)
            current_str_idx = 0
            working_variant = ""

            for (numstr, start_idx, end_idx), status in num_comb_zip:
                working_variant += self._gene_name[current_str_idx:start_idx]
                if status == "keep":
                    working_variant += f"{numstr}-1"
                else:
                    working_variant += numstr
                current_str_idx = end_idx

            working_variant += self._gene_name[current_str_idx:]

            if working_variant != original_gene_name:
                dash1_variants.append(working_variant)

        for variant in dash1_variants:
            # specifically check _valid_gene_dictionary (don't convert from gene to subgroup)
            if variant in self._valid_gene_dictionary:
                self._gene_name = variant
                return



