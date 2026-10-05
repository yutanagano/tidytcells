import itertools
import re
from abc import ABC, abstractmethod
from typing import Dict, List, Optional

from tidytcells import _utils
from tidytcells._resources import (
    VALID_HOMOSAPIENS_MH,
    HOMOSAPIENS_MH_SYNONYMS,
    VALID_HOMOSAPIENS_MH_MRO,
    HOMOSAPIENS_MH_SYNONYMS_ALLELE_MRO,
)
from tidytcells.result import HLAGene


class HlaSymbolParser:
    MUTANT_REGEX = re.compile(r"^\s*(\S+)\s+(.+?)\s+mutant\s*$", re.IGNORECASE)

    cleaned_symbol: str
    gene_name: str
    allele_designation: List[str]
    mutation: Optional[str]

    def __init__(self, hla_symbol: str) -> None:
        hla_symbol = self._parse_mutation(hla_symbol)
        self.cleaned_symbol = _utils.clean_and_uppercase(hla_symbol)
        self._parse_gene_and_allele(self.cleaned_symbol)

    def _parse_mutation(self, hla_symbol: str) -> str:
        self.mutation = None

        mutant_match = self.MUTANT_REGEX.match(hla_symbol)
        if mutant_match:
            hla_symbol, self.mutation = mutant_match.groups()

        return hla_symbol

    def _parse_gene_and_allele(self, hla_symbol: str) -> None:
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
            f"{int(d):02}" if d.isdigit() and len(d) <= 3 else d
            for d in allele_designation.split(":")
        ]


class HlaSymbolStandardizer(ABC):
    """
    Abstract base standardizer class.
    """

    @property
    @abstractmethod
    def _valid_symbols(self) -> Dict[str, Dict]:
        pass

    @property
    @abstractmethod
    def _gene_synonyms(self) -> Dict[str, str]:
        pass

    @property
    @abstractmethod
    def _allele_synonyms(self) -> Dict[str, str]:
        pass

    def __init__(self, symbol: str) -> None:
        self.original_symbol = symbol
        self._parse_hla_symbol(symbol)
        self._resolve_errors()
        self._compile_result()

    def _parse_hla_symbol(self, hla_symbol: str) -> None:
        parsed_hla_symbol = HlaSymbolParser(hla_symbol)
        self._mutation = parsed_hla_symbol.mutation

        if parsed_hla_symbol.cleaned_symbol in self._allele_synonyms:
            parsed_hla_symbol = HlaSymbolParser(
                self._allele_synonyms[parsed_hla_symbol.cleaned_symbol]
            )

        self._gene_name = parsed_hla_symbol.gene_name
        self._allele_designation = parsed_hla_symbol.allele_designation

    def _resolve_errors(self, skip_add1_section: bool = False) -> None:
        if self.get_reason_why_invalid() is None:
            return

        if self._is_synonym():
            self._gene_name = self._gene_synonyms[self._gene_name]
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

        if not skip_add1_section:
            original = self._gene_name
            if original.endswith("1"):
                self._gene_name = original[:-1]
            else:
                self._gene_name = original + "1"
            if self.get_reason_why_invalid() is None:
                return
            self._gene_name = original

    def _is_synonym(self) -> bool:
        return self._gene_name in self._gene_synonyms

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
        if not self._allele_designation or not self._allele_designation[0].isdigit():
            return

        original = self._allele_designation

        for fields in self._split_into_designator_fields(original[0]):
            self._allele_designation = fields + original[1:]
            if self.get_reason_why_invalid() is None:
                return

        self._allele_designation = original

    @classmethod
    def _split_into_designator_fields(cls, digits: str, max_fields: int = 4) -> List[List[str]]:
        candidates = [f for f in cls._all_designator_field_splits(digits) if 2 <= len(f) <= max_fields]
        return sorted(candidates, key=cls._field_lengths)

    @classmethod
    def _all_designator_field_splits(cls, digits: str) -> List[List[str]]:
        if not digits:
            return [[]]

        splits = []
        for n in (2, 3):
            field = digits[:n]
            if len(field) == n and not (n == 3 and field.startswith("0")):
                for tail in cls._all_designator_field_splits(digits[n:]):
                    splits.append([field] + tail)

        return splits

    @staticmethod
    def _field_lengths(fields: List[str]) -> List[int]:
        return [len(field) for field in fields]

    def _try_different_amounts_of_leading_zeros_in_first_2_allele_designators(
        self,
    ) -> None:
        if not all(ad.isdigit() for ad in self._allele_designation[:2]):
            return

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

    def get_reason_why_invalid(self) -> Optional[str]:
        if self._gene_name == "B2M" and not self._allele_designation:
            return None

        if not self._gene_name in self._valid_symbols:
            return "Unrecognized gene name"

        # Verify allele designators up to the level of the protein (or G/P)
        allele_designation = self._allele_designation.copy()
        if not self._is_group():
            allele_designation = allele_designation[:2]
        current_root = self._valid_symbols[self._gene_name]

        while len(allele_designation) > 0:
            try:
                current_root = current_root[allele_designation.pop(0)]
            except KeyError:
                return "Nonexistent allele for recognized gene"

        # If there are designator fields past the protein level, just make sure
        # they look like legitimate designator field values
        if not self._is_group() and len(self._allele_designation) > 2:
            further_designators = self._allele_designation[2:]

            if len(further_designators) > 2:
                return "Too many allele designators"

            for field in further_designators:
                if not field.isdigit():
                    return "Non-numerical allele designators"

                if len(field) < 2:
                    return "Non-2-digit allele designators"

        return None

    def _is_group(self) -> bool:
        if not self._allele_designation:
            return False

        return self._allele_designation[-1].endswith("G") or self._allele_designation[
            -1
        ].endswith("P")


    def _compile_result(self):
        self.result = HLAGene(original_input=self.original_symbol,
                              error=self.get_reason_why_invalid(),
                              gene_name=self._gene_name,
                              allele_designation=self._allele_designation,
                              mutation=self._mutation)


class ImgtHlaSymbolStandardizer(HlaSymbolStandardizer):
    _valid_symbols = VALID_HOMOSAPIENS_MH
    _gene_synonyms = HOMOSAPIENS_MH_SYNONYMS
    _allele_synonyms = {}


class MroHlaSymbolStandardizer(HlaSymbolStandardizer):
    _valid_symbols = VALID_HOMOSAPIENS_MH_MRO
    _gene_synonyms = {}
    _allele_synonyms = HOMOSAPIENS_MH_SYNONYMS_ALLELE_MRO