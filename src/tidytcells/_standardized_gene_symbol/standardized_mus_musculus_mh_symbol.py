import re
from abc import ABC, abstractmethod
from collections import defaultdict
from typing import Dict, Optional

from tidytcells import _utils
from tidytcells._resources import (
    VALID_MUSMUSCULUS_MH,
    MUSMUSCULUS_MH_SYNONYMS,
    VALID_MUSMUSCULUS_MH_MRO,
    MUSMUSCULUS_MH_SYNONYMS_ALLELE_MRO,
)
from tidytcells.result import MhGene


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


class MroMhSymbolParser:
    H2_PREFIX_REGEX = re.compile(r"^H-?2-?", re.IGNORECASE)

    cleaned_symbol: str
    has_h2_prefix: bool
    core: str

    def __init__(self, mh_symbol: str) -> None:
        self.cleaned_symbol = re.sub(r"\s+", "", mh_symbol)
        self.has_h2_prefix = self.H2_PREFIX_REGEX.match(self.cleaned_symbol) is not None
        self.core = self.strip_h2_prefix_and_hyphens(self.cleaned_symbol)

    @classmethod
    def strip_h2_prefix_and_hyphens(cls, mh_symbol: str) -> str:
        return cls.H2_PREFIX_REGEX.sub("", mh_symbol).replace("-", "")


class MusMusculusMhSymbolStandardizer(ABC):
    """
    Abstract base standardizer class.
    """

    @property
    @abstractmethod
    def _valid_symbols(self) -> Dict[str, Dict]:
        pass

    @property
    @abstractmethod
    def _synonyms(self) -> Dict[str, str]:
        pass

    def __init__(self, symbol: str) -> None:
        self.original_symbol = symbol
        self._allele_designation = None
        self._standardize(symbol)
        self._compile_result()

    @abstractmethod
    def _standardize(self, symbol: str) -> None:
        pass

    def get_reason_why_invalid(self, enforce_functional: bool = False) -> Optional[str]:
        if not self._gene_name in self._valid_symbols:
            return "Unrecognized gene name"

        return None

    def _compile_result(self):
        self.result = MhGene(original_input=self.original_symbol,
                             error=self.get_reason_why_invalid(),
                             gene_name=self._gene_name,
                             allele_designation=self._allele_designation,
                             species="musmusculus")


class ImgtMusMusculusMhSymbolStandardizer(MusMusculusMhSymbolStandardizer):
    _valid_symbols = VALID_MUSMUSCULUS_MH
    _synonyms = MUSMUSCULUS_MH_SYNONYMS

    def _standardize(self, symbol: str) -> None:
        self._parse_mh_symbol(symbol)
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
            self._gene_name = self._synonyms[self._gene_name.replace("-", "")]
            if self.get_reason_why_invalid() is None:
                return

    def _is_synonym(self) -> bool:
        return self._gene_name.replace("-", "") in self._synonyms

    def compile(self, precision: str = "allele") -> str:
        if precision == "allele" and self._allele_designation:
            return f"{self._gene_name}*{self._allele_designation}"

        return self._gene_name


def build_prefixed_index(valid_symbols: Dict, synonyms: Dict[str, str]) -> Dict[str, str]:
    index = {MroMhSymbolParser.strip_h2_prefix_and_hyphens(s).upper(): s for s in valid_symbols}

    synonym_targets = defaultdict(set)
    for synonym, target in synonyms.items():
        synonym_targets[MroMhSymbolParser.strip_h2_prefix_and_hyphens(synonym).upper()].add(target)

    for core, targets in synonym_targets.items():
        if core not in index and len(targets) == 1:
            index[core] = targets.pop()

    return index


class MroMusMusculusMhSymbolStandardizer(MusMusculusMhSymbolStandardizer):
    _valid_symbols = VALID_MUSMUSCULUS_MH_MRO
    _synonyms = MUSMUSCULUS_MH_SYNONYMS_ALLELE_MRO

    _exact_index = {**{s.upper(): s for s in _valid_symbols}, **_synonyms}
    _prefixed_index = build_prefixed_index(_valid_symbols, _synonyms)
    _bare_index = {MroMhSymbolParser.strip_h2_prefix_and_hyphens(s): s for s in _valid_symbols}

    def _standardize(self, symbol: str) -> None:
        parsed_mh_symbol = MroMhSymbolParser(symbol)
        self._gene_name = self._resolve(parsed_mh_symbol)

    def _resolve(self, parsed_mh_symbol: MroMhSymbolParser) -> str:
        cleaned_symbol = parsed_mh_symbol.cleaned_symbol

        if cleaned_symbol.upper() in self._exact_index:
            return self._exact_index[cleaned_symbol.upper()]

        # without an H2 prefix, case is needed to tell haplotypes apart from other names (e.g. "Dr" vs "DR")
        if parsed_mh_symbol.has_h2_prefix:
            return self._prefixed_index.get(parsed_mh_symbol.core.upper(), cleaned_symbol)

        return self._bare_index.get(parsed_mh_symbol.core, cleaned_symbol)