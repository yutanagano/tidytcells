import re
from collections import defaultdict
from typing import Dict, Optional

from tidytcells._standardized_gene_symbol.mutation import interpret_remainder
from tidytcells.result import MhGene


def strip_prefix_and_hyphens(symbol: str, prefix_regex: re.Pattern) -> str:
    return prefix_regex.sub("", symbol).replace("-", "")


class MhHaplotypeSymbolParser:
    cleaned_symbol: str
    token_symbol: str
    remainder: str

    def __init__(self, symbol: str) -> None:
        symbol = symbol.strip()
        parts = symbol.split(None, 1)
        self.token_symbol = parts[0]
        self.remainder = parts[1].strip() if len(parts) > 1 else ""
        self.cleaned_symbol = re.sub(r"\s+", "", symbol)


def build_prefixed_index(
    valid_symbols: Dict, allele_synonyms: Dict[str, str], prefix_regex: re.Pattern
) -> Dict[str, str]:
    index = {
        strip_prefix_and_hyphens(s, prefix_regex).upper(): s for s in valid_symbols
    }

    synonym_targets = defaultdict(set)
    for synonym, target in allele_synonyms.items():
        synonym_targets[strip_prefix_and_hyphens(synonym, prefix_regex).upper()].add(
            target
        )

    for core, targets in synonym_targets.items():
        if core not in index and len(targets) == 1:
            index[core] = targets.pop()

    return index


class MhHaplotypeSymbolStandardizer:
    """
    Standardizer for MH symbols named by locus and haplotype (e.g. H2-Kb, H2-IAb).
    """

    def __init__(
        self,
        symbol: str,
        species: str,
        mh_database: str,
        valid_symbols: Dict[str, Dict],
        allele_synonyms: Dict[str, str],
        gene_synonyms=None,
    ):
        self.original_symbol = symbol
        self.species = species
        self.mh_database = mh_database
        self._valid_symbols = valid_symbols
        self._allele_synonyms = allele_synonyms

        if species == "musmusculus":
            self._prefix_regex = re.compile(r"^H-?2-?", re.IGNORECASE)
        elif species == "rattusnorvegicus":
            self._prefix_regex = re.compile(r"^RT-?1-?", re.IGNORECASE)

        self._exact_index = {**allele_synonyms, **{s.upper(): s for s in valid_symbols}}
        self._prefixed_index = build_prefixed_index(
            valid_symbols, allele_synonyms, self._prefix_regex
        )
        self._bare_index = {
            strip_prefix_and_hyphens(s, self._prefix_regex): s for s in valid_symbols
        }

        parsed_symbol = MhHaplotypeSymbolParser(symbol)

        # a valid haplotype symbol never contains whitespace, so if the whole
        # input (whitespace removed) resolves, the whitespace was a typo
        # (e.g. "H-2 IEk" -> "H2-IEk") and there is no remainder to interpret
        whole_input = self._resolve(parsed_symbol.cleaned_symbol)
        if whole_input in self._valid_symbols:
            self._gene_name = whole_input
            remainder = ""
        else:
            self._gene_name = self._resolve(parsed_symbol.token_symbol)
            remainder = parsed_symbol.remainder.upper()

        self._mutation, self._mutation_error = interpret_remainder(remainder)
        self._trailing_text = (
            remainder
            if remainder and self._mutation is None and self._mutation_error is None
            else None
        )

        self._compile_result()

    def _resolve(self, cleaned_symbol: str) -> str:
        if cleaned_symbol.upper() in self._exact_index:
            return self._exact_index[cleaned_symbol.upper()]

        core = strip_prefix_and_hyphens(cleaned_symbol, self._prefix_regex)

        # without a prefix, case is needed to tell haplotypes apart from other names (e.g. "Dr" vs "DR")
        if self._prefix_regex.match(cleaned_symbol) is not None:
            return self._prefixed_index.get(core.upper(), cleaned_symbol)

        return self._bare_index.get(core, cleaned_symbol)

    def get_reason_why_invalid(self) -> Optional[str]:
        if self._gene_name not in self._valid_symbols:
            return "Unrecognized gene name"

        # the remainder after the symbol, if any, must be a valid mutation
        # specification
        if self._mutation_error is not None:
            return self._mutation_error

        if self._trailing_text is not None:
            return f"Unexpected trailing text: {self._trailing_text}"

        return None

    def _compile_result(self) -> None:
        self.result = MhGene(
            original_input=self.original_symbol,
            error=self.get_reason_why_invalid(),
            gene_name=self._gene_name,
            species=self.species,
            mutation=self._mutation,
            mh_database=self.mh_database,
        )
