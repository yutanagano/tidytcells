import re
import warnings

from tidytcells._resources import MH_SYMBOLS_BY_SP_DB


def query_symbols(
    precision: str,
    contains_pattern: str | None = None,
    *,
    species: str,
    mh_database: str,
    warning: str | None = None,
) -> frozenset[str]:
    """
    Query all known MH gene/allele symbols for a species and MH database at
    the requested precision, optionally filtered by a regular expression.

    Genes with no allele-level data contribute their bare gene name, since
    that is the finest resolution the data supports.
    """
    valid_symbols = _valid_symbols(species, mh_database)

    if precision == "allele" and warning is not None:
        warnings.warn(warning)

    pattern = re.compile(contains_pattern) if contains_pattern is not None else None

    def keep(symbol: str) -> bool:
        return pattern is None or pattern.search(symbol) is not None

    query_results = []

    for gene_symbol, allele_dictionary in valid_symbols.items():
        if precision == "gene" or not allele_dictionary:
            if keep(gene_symbol):
                query_results.append(gene_symbol)
            continue

        for first_allele_designation in allele_dictionary:
            for second_allele_designation in allele_dictionary[
                first_allele_designation
            ]:
                if not second_allele_designation.isdigit():
                    continue

                symbol = (
                    gene_symbol
                    + "*"
                    + first_allele_designation
                    + ":"
                    + second_allele_designation
                )
                if keep(symbol):
                    query_results.append(symbol)

    return frozenset(query_results)


def _valid_symbols(species: str, mh_database: str) -> dict:
    """
    The valid-symbol tree for a species and MH database.
    """
    return MH_SYMBOLS_BY_SP_DB[species][mh_database]["symbols"]
