import re
from abc import ABC, abstractmethod
from collections import defaultdict
from typing import Dict, Optional

from tidytcells._resources import *
from tidytcells.result import MhGene


def strip_prefix_and_hyphens(symbol: str, prefix_regex: re.Pattern) -> str:
    return prefix_regex.sub("", symbol).replace("-", "")


class MhHaplotypeSymbolParser:
    MUTANT_REGEX = re.compile(r"^\s*(\S+)\s+(.+?)\s+mutant\s*$", re.IGNORECASE)

    cleaned_symbol: str
    has_prefix: bool
    core: str
    mutation: Optional[str]

    def __init__(self, symbol: str, prefix_regex: re.Pattern) -> None:
        symbol = self._parse_mutation(symbol)
        self.cleaned_symbol = re.sub(r"\s+", "", symbol)
        self.has_prefix = prefix_regex.match(self.cleaned_symbol) is not None
        self.core = strip_prefix_and_hyphens(self.cleaned_symbol, prefix_regex)

    def _parse_mutation(self, symbol: str) -> str:
        self.mutation = None

        mutant_match = self.MUTANT_REGEX.match(symbol)
        if mutant_match:
            symbol, self.mutation = mutant_match.groups()

        return symbol


def build_prefixed_index(valid_symbols: Dict, allele_synonyms: Dict[str, str], prefix_regex: re.Pattern) -> Dict[str, str]:
    index = {strip_prefix_and_hyphens(s, prefix_regex).upper(): s for s in valid_symbols}

    synonym_targets = defaultdict(set)
    for synonym, target in allele_synonyms.items():
        synonym_targets[strip_prefix_and_hyphens(synonym, prefix_regex).upper()].add(target)

    for core, targets in synonym_targets.items():
        if core not in index and len(targets) == 1:
            index[core] = targets.pop()

    return index


class MhHaplotypeSymbolStandardizer:
    """
    Standardizer for MH symbols named by locus and haplotype (e.g. H2-Kb, H2-IAb).
    """

    def __init__(self, symbol: str, species: str, valid_symbols: Dict[str, Dict],
                 allele_synonyms: Dict[str, str], gene_synonyms = None):
        self.original_symbol = symbol
        self._species = species
        self._valid_symbols = valid_symbols
        self._allele_synonyms = allele_synonyms

        if species == "musmusculus":
            self._prefix_regex = re.compile(r"^H-?2-?", re.IGNORECASE)
        elif species == "rattusnorvegicus":
            self._prefix_regex = re.compile(r"^RT-?1-?", re.IGNORECASE)

        self._exact_index = {**allele_synonyms, **{s.upper(): s for s in valid_symbols}}
        self._prefixed_index = build_prefixed_index(valid_symbols, allele_synonyms, self._prefix_regex)
        self._bare_index = {strip_prefix_and_hyphens(s, self._prefix_regex): s for s in valid_symbols}

        parsed_symbol = MhHaplotypeSymbolParser(symbol, self._prefix_regex)
        self._mutation = parsed_symbol.mutation
        self._gene_name = self._resolve(parsed_symbol)

        self._compile_result()


    def _resolve(self, parsed_symbol: MhHaplotypeSymbolParser) -> str:
        cleaned_symbol = parsed_symbol.cleaned_symbol

        if cleaned_symbol.upper() in self._exact_index:
            return self._exact_index[cleaned_symbol.upper()]

        # without a prefix, case is needed to tell haplotypes apart from other names (e.g. "Dr" vs "DR")
        if parsed_symbol.has_prefix:
            return self._prefixed_index.get(parsed_symbol.core.upper(), cleaned_symbol)

        return self._bare_index.get(parsed_symbol.core, cleaned_symbol)

    def get_reason_why_invalid(self) -> Optional[str]:
        if self._gene_name not in self._valid_symbols:
            return "Unrecognized gene name"

        return None

    def _compile_result(self) -> None:
        self.result = MhGene(original_input=self.original_symbol,
                             error=self.get_reason_why_invalid(),
                             gene_name=self._gene_name,
                             species=self._species,
                             mutation=self._mutation)

