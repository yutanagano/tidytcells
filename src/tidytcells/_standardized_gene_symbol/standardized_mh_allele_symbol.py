import itertools
import re
from abc import ABC, abstractmethod
from typing import Dict, List, Optional

from tidytcells import _utils
from tidytcells.result import HLAGene, MhGene

from tidytcells._resources import *

def uppercase_tree(tree: Dict[str, Dict]) -> Dict[str, Dict]:
    return {key.upper(): uppercase_tree(subtree) for key, subtree in tree.items()}


class MhAlleleSymbolParser:
    MUTANT_REGEX = re.compile(r"^\s*(\S+)\s+(.+?)\s+mutant\s*$", re.IGNORECASE)
    SYMBOL_REGEX = re.compile(r"^([A-Z0-9\-\.\:\/]+)(\*([A-Z\d:]+))?")

    cleaned_symbol: str
    gene_name: str
    allele_designation: List[str]
    mutation: Optional[str]

    def __init__(self, symbol: str) -> None:
        symbol = self._parse_mutation(symbol)
        self.cleaned_symbol = _utils.clean_and_uppercase(symbol)
        self._parse_gene_and_allele(self._replace_periods_between_digits_with_colon(self.cleaned_symbol))

    def _parse_mutation(self, symbol: str) -> str:
        self.mutation = None

        mutant_match = self.MUTANT_REGEX.match(symbol)
        if mutant_match:
            symbol, self.mutation = mutant_match.groups()

        return symbol

    def _parse_gene_and_allele(self, symbol: str) -> None:
        parse_attempt = self.SYMBOL_REGEX.match(symbol)

        if parse_attempt:
            self.gene_name = parse_attempt.group(1)
            self.allele_designation = self._listify_allele_designation(parse_attempt.group(3))
            return

        self.gene_name = symbol
        self.allele_designation = []

    def _replace_periods_between_digits_with_colon(self, string: str) -> str:
        return re.sub(r"(?<=\d)\.(?=\d)", ":", string)

    def _listify_allele_designation(self, allele_designation: Optional[str]) -> List[str]:
        if allele_designation is None:
            return []

        allele_designation = re.sub(r"(?<=\d)[LSCAQN]$", "", allele_designation)

        return [d.zfill(2) if d.isdigit() else d for d in allele_designation.split(":")]


class MhAlleleSymbolStandardizer(ABC):
    """
    Abstract base standardizer class for MH symbols with numbered alleles (e.g. HLA-A*02:01, Mamu-A1*001:01).
    Input is matched case-insensitively; results are returned in the case used by the database.
    """

    _parser_cls = MhAlleleSymbolParser
    _gene_synonyms: Dict[str, str] = {}
    _allow_zero_padded_three_digit_fields = True

    @property
    @abstractmethod
    def _species(self) -> str:
        pass

    @property
    @abstractmethod
    def _valid_symbols(self) -> Dict[str, Dict]:
        pass

    @property
    @abstractmethod
    def _allele_synonyms(self) -> Dict[str, str]:
        pass

    def __init_subclass__(cls, **kwargs) -> None:
        super().__init_subclass__(**kwargs)

        if isinstance(cls._valid_symbols, dict):
            cls._valid_symbols_upper = uppercase_tree(cls._valid_symbols)
            cls._canonical_gene_names = {gene.upper(): gene for gene in cls._valid_symbols}
            cls._prefixes = sorted({gene.split("-")[0].upper() for gene in cls._valid_symbols if "-" in gene})

    def __init__(self, symbol: str) -> None:
        self.original_symbol = symbol
        self._parse_symbol(symbol)
        self._resolve_errors()
        self._compile_result()

    def _parse_symbol(self, symbol: str) -> None:
        parsed_symbol = self._parser_cls(symbol)
        self._mutation = parsed_symbol.mutation
        self._cleaned_symbol = parsed_symbol.cleaned_symbol
        self._gene_name = parsed_symbol.gene_name
        self._allele_designation = parsed_symbol.allele_designation

    def _resolve_errors(self) -> None:
        if self.get_reason_why_invalid() is None:
            return

        if self._cleaned_symbol in self._allele_synonyms:
            parsed_synonym = self._parser_cls(self._allele_synonyms[self._cleaned_symbol])
            self._gene_name = parsed_synonym.gene_name
            self._allele_designation = parsed_synonym.allele_designation
            if self.get_reason_why_invalid() is None:
                return

        if self._gene_name in self._gene_synonyms:
            self._gene_name = self._gene_synonyms[self._gene_name].upper()
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
        if self.get_reason_why_invalid() is None:
            return

        self._try_adding_or_removing_trailing_1_in_gene_name()

    def _resolve_common_errors(self) -> None:
        self._add_missing_prefix()

    def _add_missing_prefix(self) -> None:
        if any(self._gene_name.startswith(f"{prefix}-") for prefix in self._prefixes):
            return

        original = self._gene_name

        for prefix in self._prefixes:
            self._gene_name = f"{prefix}-{original}"
            if self._gene_name in self._valid_symbols_upper:
                return

        self._gene_name = f"{self._prefixes[0]}-{original}" if len(self._prefixes) == 1 else original

    def _handle_forgotten_asterisk(self) -> None:
        if self._allele_designation:
            return

        for split_idx in range(len(self._gene_name) - 1, 0, -1):
            gene, designation = self._gene_name[:split_idx], self._gene_name[split_idx:]

            if gene in self._valid_symbols_upper and re.fullmatch(r"[\d:]+", designation):
                self._gene_name = gene
                self._allele_designation = designation.split(":")
                return

    def _handle_forgotten_colon_between_first_and_second_allele_designator(self) -> None:
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
            zero_padded = n == 3 and field.startswith("0")

            if len(field) == n and (cls._allow_zero_padded_three_digit_fields or not zero_padded):
                for tail in cls._all_designator_field_splits(digits[n:]):
                    splits.append([field] + tail)

        return splits

    @staticmethod
    def _field_lengths(fields: List[str]) -> List[int]:
        return [len(field) for field in fields]

    def _try_different_amounts_of_leading_zeros_in_first_2_allele_designators(self) -> None:
        original = self._allele_designation
        reformatted_allele_designators = [
            [f"{int(ad):02}", f"{int(ad):03}"] if ad.isdigit() else [ad] for ad in original[:2]
        ]

        for new_designators in itertools.product(*reformatted_allele_designators):
            self._allele_designation = list(new_designators) + original[2:]
            if self.get_reason_why_invalid() is None:
                return

        self._allele_designation = original

    def _try_adding_or_removing_trailing_1_in_gene_name(self) -> None:
        original = self._gene_name
        self._gene_name = original[:-1] if original.endswith("1") else original + "1"

        if self.get_reason_why_invalid() is not None:
            self._gene_name = original

    def get_reason_why_invalid(self) -> Optional[str]:
        if self._gene_name not in self._valid_symbols_upper:
            return "Unrecognized gene name"

        # Verify allele designators up to the level of the protein (or G/P)
        allele_designation = self._allele_designation.copy()
        if not self._is_group():
            allele_designation = allele_designation[:2]
        current_root = self._valid_symbols_upper[self._gene_name]

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

        return self._allele_designation[-1].endswith(("G", "P"))

    def _canonical_gene_name(self) -> str:
        return self._canonical_gene_names.get(self._gene_name, self._gene_name)

    def _canonical_allele_designation(self) -> List[str]:
        node = self._valid_symbols.get(self._canonical_gene_name(), dict())
        canonical = []

        for field in self._allele_designation:
            canonical_field = next((key for key in node if key.upper() == field), field)
            canonical.append(canonical_field)
            node = node.get(canonical_field, dict())

        return canonical

    def _compile_result(self) -> None:
        self.result = MhGene(original_input=self.original_symbol,
                             error=self.get_reason_why_invalid(),
                             gene_name=self._canonical_gene_name(),
                             allele_designation=self._canonical_allele_designation(),
                             species=self._species,
                             mutation=self._mutation)


class HlaSymbolParser(MhAlleleSymbolParser):
    SYMBOL_REGEX = re.compile(r"^([A-Z0-9\-\.\:\/]+)(\*([\d:]+G?P?)[LSCAQN]?)?")
    HLA_CLASS_II_REGEX = re.compile(r"^((HLA-)?(D[PQ][AB]|DRB|TAP)\d)(\*?([\d:]+G?P?)[LSCAQN]?)?")

    def _parse_gene_and_allele(self, symbol: str) -> None:
        if symbol == "B2M":
            self.gene_name = "B2M"
            self.allele_designation = []
            return

        parse_attempt = self.HLA_CLASS_II_REGEX.match(symbol)
        if parse_attempt:
            self.gene_name = parse_attempt.group(1)
            self.allele_designation = self._listify_allele_designation(parse_attempt.group(5))
            return

        super()._parse_gene_and_allele(symbol)


class HlaSymbolStandardizer(MhAlleleSymbolStandardizer):
    """
    Abstract base standardizer class for HLA symbols.
    """

    _species = "homosapiens"
    _parser_cls = HlaSymbolParser
    _allow_zero_padded_three_digit_fields = False

    def _resolve_common_errors(self) -> None:
        super()._resolve_common_errors()
        self._gene_name = self._gene_name.replace("CW", "C")

    def get_reason_why_invalid(self) -> Optional[str]:
        if self._gene_name == "B2M" and not self._allele_designation:
            return None

        return super().get_reason_why_invalid()

    def _compile_result(self) -> None:
        self.result = HLAGene(original_input=self.original_symbol,
                              error=self.get_reason_why_invalid(),
                              gene_name=self._canonical_gene_name(),
                              allele_designation=self._canonical_allele_designation(),
                              mutation=self._mutation)


class ImgtHlaSymbolStandardizer(HlaSymbolStandardizer):
    _valid_symbols = HOMOSAPIENS_MH_VALID_IPD
    _gene_synonyms = HOMOSAPIENS_MH_SYNONYMS_GENE_HGNC_IPD
    _allele_synonyms = {}

class MroHlaSymbolStandardizer(HlaSymbolStandardizer):
    _valid_symbols = HOMOSAPIENS_MH_VALID_PROTEIN_MRO
    _gene_synonyms = HOMOSAPIENS_MH_SYNONYMS_GENE_HGNC_MRO
    _allele_synonyms = HOMOSAPIENS_MH_SYNONYMS_ALLELE_MRO



class MroBoLaSymbolStandardizer(MhAlleleSymbolStandardizer):
    _species = "bostaurus"
    _parser_cls = MhAlleleSymbolParser
    _allow_zero_padded_three_digit_fields = True
    _valid_symbols = BOSTAURUS_MH_VALID_PROTEIN_MRO
    _gene_synonyms = {}
    _allele_synonyms = BOSTAURUS_MH_SYNONYMS_ALLELE_MRO

class MroMamuSymbolStandardizer(MhAlleleSymbolStandardizer):
    _species = "macacamulatta"
    _parser_cls = MhAlleleSymbolParser
    _allow_zero_padded_three_digit_fields = True
    _valid_symbols = MACACAMULATTA_MH_VALID_PROTEIN_MRO
    _gene_synonyms = {}
    _allele_synonyms = MACACAMULATTA_MH_SYNONYMS_ALLELE_MRO

class MroPatrSymbolStandardizer(MhAlleleSymbolStandardizer):
    _species = "pantroglodytes"
    _parser_cls = MhAlleleSymbolParser
    _allow_zero_padded_three_digit_fields = False
    _valid_symbols = PANTROGLODYTES_MH_VALID_PROTEIN_MRO
    _gene_synonyms = {}
    _allele_synonyms = PANTROGLODYTES_MH_SYNONYMS_ALLELE_MRO

class MroSlaSymbolStandardizer(MhAlleleSymbolStandardizer):
    _species = "susscrofa"
    _parser_cls = MhAlleleSymbolParser
    _allow_zero_padded_three_digit_fields = False
    _valid_symbols = SUSSCROFA_MH_VALID_PROTEIN_MRO
    _gene_synonyms = {}
    _allele_synonyms = SUSSCROFA_MH_SYNONYMS_ALLELE_MRO




