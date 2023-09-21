import itertools
import re
from typing import List, Optional

from tidytcells import _utils
from tidytcells._standardized_gene_symbol import StandardizedGeneSymbol
from tidytcells._resources import VALID_HOMOSAPIENS_MH, HOMOSAPIENS_MH_SYNONYMS


class HlaSymbolParser:
    gene_name: str
    allele_designation: List[str]

    def __init__(self, hla_symbol: str) -> None:
        if hla_symbol == "B2M":
            self.gene_name = "B2M"
            self.allele_designation = []
            return

        hla_symbol = self._replace_periods_between_digits_with_colon(hla_symbol)

        parse_attempt_1 = re.match(
            r"^((HLA-)?(D[PQ][AB]|DRB|TAP)\d)(\*?([\d:]+G?P?)[LSCAQN]?)?", hla_symbol
        )
        if parse_attempt_1:
            self.gene_name = parse_attempt_1.group(1)
            self.allele_designation = self._listify_allele_designation(
                parse_attempt_1.group(5)
            )
            return

        parse_attempt_2 = re.match(
            r"^([A-Z0-9\-\.\:\/]+)(\*([\d:]+G?P?)[LSCAQN]?)?", hla_symbol
        )
        if parse_attempt_2:
            self.gene_name = parse_attempt_2.group(1)
            self.allele_designation = self._listify_allele_designation(
                parse_attempt_2.group(3)
            )
            return

        self.gene_name = hla_symbol
        self.allele_designation = []

    def _replace_periods_between_digits_with_colon(self, string: str) -> str:
        return re.sub(r"(?<=\d)\.(?=\d)", ":", string)

    def _listify_allele_designation(
        self, allele_designation: Optional[str]
    ) -> List[str]:
        if allele_designation is None:
            return []

        return [
            f"{int(d):02}" if d.isdigit() else d for d in allele_designation.split(":")
        ]


class StandardizedHlaSymbol(StandardizedGeneSymbol):
    def __init__(self, gene_symbol: str) -> None:
        self._parse_hla_symbol(gene_symbol)
        self._resolve_errors()

    def _parse_hla_symbol(self, hla_symbol: str) -> None:
        cleaned_hla_symbol = _utils.clean_and_uppercase(hla_symbol)
        parsed_hla_symbol = HlaSymbolParser(cleaned_hla_symbol)
        self._gene_name = parsed_hla_symbol.gene_name
        self._allele_designation = parsed_hla_symbol.allele_designation

    def _resolve_errors(self) -> None:
        if self.get_reason_why_invalid() is None:
            return

        if self._is_synonym():
            self._gene_name = HOMOSAPIENS_MH_SYNONYMS[self._gene_name]
            if self.get_reason_why_invalid() is None:
                return

        self._resolve_common_errors()
        if self.get_reason_why_invalid() is None:
            return

        self._handle_forgotten_asterisk()
        if self.get_reason_why_invalid() is None:
            return

        self._handle_forgotten_colon_between_first_and_second_allele_designator()
        if self.get_reason_why_invalid() is None:
            return

        self._try_different_amounts_of_leading_zeros_in_first_2_allele_designators()

    def _is_synonym(self) -> bool:
        return self._gene_name in HOMOSAPIENS_MH_SYNONYMS

    def _resolve_common_errors(self) -> None:
        if not self._gene_name.startswith("HLA-"):
            self._gene_name = "HLA-" + self._gene_name
        self._gene_name = self._gene_name.replace("CW", "C")

    def _handle_forgotten_asterisk(self) -> None:
        if not self._allele_designation:
            m = re.match(r"^(HLA-[A-Z]+)([\d:]+G?P?)$", self._gene_name)
            if m:
                self._gene_name = m.group(1)
                self._allele_designation = m.group(2).split(":")

    def _handle_forgotten_colon_between_first_and_second_allele_designator(
        self,
    ) -> None:
        if self._allele_designation and len(self._allele_designation[0]) == 4:
            self._allele_designation = [
                self._allele_designation[0][:2],
                self._allele_designation[0][2:],
            ] + self._allele_designation[1:]

    def _try_different_amounts_of_leading_zeros_in_first_2_allele_designators(
        self,
    ) -> None:
        original = self._allele_designation
        first_two_allele_designators = [int(ad) for ad in self._allele_designation[:2]]
        reformatted_allele_designators = [
            [f"{ad:02}", f"{ad:03}"] for ad in first_two_allele_designators
        ]

        for new_designators in itertools.product(*reformatted_allele_designators):
            self._allele_designation = (
                list(new_designators) + self._allele_designation[2:]
            )
            if self.get_reason_why_invalid() is None:
                return

        self._allele_designation = original

    def get_reason_why_invalid(self, enforce_functional: bool = False) -> Optional[str]:
        if self._gene_name == "B2M" and not self._allele_designation:
            return None

        if not self._gene_name in VALID_HOMOSAPIENS_MH:
            return "unrecognised gene name"

        # Verify allele designators up to the level of the protein (or G/P)
        allele_designation = self._allele_designation.copy()
        if not self._is_group():
            allele_designation = allele_designation[:2]
        current_root = VALID_HOMOSAPIENS_MH[self._gene_name]

        while len(allele_designation) > 0:
            try:
                current_root = current_root[allele_designation.pop(0)]
            except KeyError:
                return "nonexistent allele for recognised gene"

        # If there are designator fields past the protein level, just make sure
        # they look like legitimate designator field values
        if not self._is_group() and len(self._allele_designation) > 2:
            further_designators = self._allele_designation[2:]

            if len(further_designators) > 2:
                return "too many allele designators"

            for field in further_designators:
                if not field.isdigit():
                    return "non-numerical allele designators"

                if len(field) < 2:
                    return "non-2-digit allele designators"

        return None

    def _is_group(self) -> bool:
        if not self._allele_designation:
            return False

        return self._allele_designation[-1].endswith("G") or self._allele_designation[
            -1
        ].endswith("P")

    def compile(self, precision: str = "allele") -> str:
        if self._allele_designation:
            if precision == "allele":
                return f'{self._gene_name}*{":".join(self._allele_designation)}'

            if precision == "protein":
                return f'{self._gene_name}*{":".join(self._allele_designation[:2])}'

        return self._gene_name
