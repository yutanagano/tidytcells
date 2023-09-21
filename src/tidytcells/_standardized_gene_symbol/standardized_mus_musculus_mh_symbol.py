import re
from typing import Optional

from tidytcells import _utils
from tidytcells._resources import VALID_MUSMUSCULUS_MH, MUSMUSCULUS_MH_SYNONYMS
from tidytcells._standardized_gene_symbol import StandardizedGeneSymbol


class MhSymbolParser:
    gene_name: str
    allele_designation: str

    def __init__(self, mh_symbol: str) -> None:
        parse_attempt = re.match(r"^([A-Z0-9\-\.\(\)\/]+)(\*(\d+))?", mh_symbol)

        if parse_attempt:
            self.gene_name = parse_attempt.group(1)
            self.allele_designation = (
                None
                if parse_attempt.group(3) is None
                else f"{int(parse_attempt.group(3)):02}"
            )
        else:
            self.gene_name = mh_symbol
            self.allele_designation = None


class StandardizedMusMusculusMhSymbol(StandardizedGeneSymbol):
    def __init__(self, gene_symbol: str) -> None:
        self._parse_mh_symbol(gene_symbol)
        self._resolve_errors()

    def _parse_mh_symbol(self, mh_symbol: str) -> None:
        cleaned_mh_symbol = _utils.clean_and_uppercase(mh_symbol)
        parsed_mh_symbol = MhSymbolParser(cleaned_mh_symbol)
        self._gene_name = parsed_mh_symbol.gene_name
        self._allele_designation = parsed_mh_symbol.allele_designation

    def _resolve_errors(self) -> None:
        if self.get_reason_why_invalid() is None:
            return

        if self._is_synonym():
            self._gene_name = MUSMUSCULUS_MH_SYNONYMS[self._gene_name.replace("-", "")]
            if self.get_reason_why_invalid() is None:
                return

    def _is_synonym(self) -> bool:
        return self._gene_name.replace("-", "") in MUSMUSCULUS_MH_SYNONYMS

    def get_reason_why_invalid(self, enforce_functional: bool = False) -> Optional[str]:
        if not self._gene_name in VALID_MUSMUSCULUS_MH:
            return "unrecognised gene name"

        return None

    def compile(self, precision: str = "allele") -> str:
        if precision == "allele" and self._allele_designation:
            return f"{self._gene_name}*{self._allele_designation}"

        return self._gene_name
