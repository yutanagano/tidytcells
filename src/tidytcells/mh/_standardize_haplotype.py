import re
from collections import defaultdict

from tidytcells._resources import MH_SYMBOLS_BY_SP_DB
from tidytcells._standardized_gene_symbol.mutation import interpret_remainder
from tidytcells.result import MhGene


def standardize_haplotype(
    symbol: str,
    *,
    species: str,
    mh_database: str,
    prefix_regex: re.Pattern,
) -> MhGene:
    """
    Standardize an MH symbol named by locus and haplotype (e.g. H2-Kb, H2-IAb).
    """
    valid_symbols, allele_synonyms = _mh_data(species, mh_database)
    exact_index = {**allele_synonyms, **{s.upper(): s for s in valid_symbols}}
    prefixed_index = build_prefixed_index(valid_symbols, allele_synonyms, prefix_regex)
    bare_index = {strip_prefix_and_hyphens(s, prefix_regex): s for s in valid_symbols}

    stripped_symbol = symbol.strip()
    parts = stripped_symbol.split(None, 1)
    token_symbol = parts[0]
    remainder = parts[1].strip() if len(parts) > 1 else ""
    cleaned_symbol = re.sub(r"\s+", "", stripped_symbol)

    # a valid haplotype symbol never contains whitespace, so if the whole
    # input (whitespace removed) resolves, the whitespace was a typo
    # (e.g. "H-2 IEk" -> "H2-IEk") and there is no remainder to interpret
    whole_input = _resolve(
        cleaned_symbol, exact_index, prefixed_index, bare_index, prefix_regex
    )
    if whole_input in valid_symbols:
        gene_name = whole_input
        remainder = ""
    else:
        gene_name = _resolve(
            token_symbol, exact_index, prefixed_index, bare_index, prefix_regex
        )
        remainder = remainder.upper()

    # only if the symbol part resolves to a standard symbol is the
    # trailing text examined as a mutation specification
    mutation = None
    if gene_name not in valid_symbols:
        error: str | None = "Unrecognized gene name"
    else:
        mutation, mutation_error = interpret_remainder(remainder)
        if mutation_error is not None:
            error = mutation_error
        elif mutation is None and remainder:
            error = f"Unexpected trailing text: {remainder}"
        else:
            error = None

    return MhGene(
        original_input=symbol,
        error=error,
        gene_name=gene_name,
        species=species,
        mutation=mutation,
        mh_database=mh_database,
    )


def _mh_data(species: str, mh_database: str) -> tuple[dict[str, dict], dict[str, str]]:
    """
    The valid-symbol and synonym data for a species and MH database.
    """
    data = MH_SYMBOLS_BY_SP_DB[species][mh_database]
    return data["symbols"], data.get("synonyms_allele", {})


def _resolve(
    cleaned_symbol: str,
    exact_index: dict[str, str],
    prefixed_index: dict[str, str],
    bare_index: dict[str, str],
    prefix_regex: re.Pattern,
) -> str:
    if cleaned_symbol.upper() in exact_index:
        return exact_index[cleaned_symbol.upper()]

    core = strip_prefix_and_hyphens(cleaned_symbol, prefix_regex)

    # without a prefix, case is needed to tell haplotypes apart from other names (e.g. "Dr" vs "DR")
    if prefix_regex.match(cleaned_symbol) is not None:
        return prefixed_index.get(core.upper(), cleaned_symbol)

    return bare_index.get(core, cleaned_symbol)


def build_prefixed_index(
    valid_symbols: dict, allele_synonyms: dict[str, str], prefix_regex: re.Pattern
) -> dict[str, str]:
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


def strip_prefix_and_hyphens(symbol: str, prefix_regex: re.Pattern) -> str:
    return prefix_regex.sub("", symbol).replace("-", "")
