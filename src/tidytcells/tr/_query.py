import re
from typing import Dict, FrozenSet, Optional, Type

from tidytcells import _utils
from tidytcells._utils import Parameter
from tidytcells._query_engine import (
    QueryEngine,
    HomoSapiensTrQueryEngine,
    MusMusculusTrQueryEngine,
)


QUERY_ENGINES: Dict[str, Type[QueryEngine]] = {
    "homosapiens": HomoSapiensTrQueryEngine,
    "musmusculus": MusMusculusTrQueryEngine,
}


def query(
    species: str = "homosapiens",
    precision: str = "allele",
    functionality: str = "any",
    contains_pattern: Optional[str] = None,
) -> FrozenSet[str]:
    """
    Query the list of all known TR genes/alleles.

    .. topic:: Supported species

        - ``"homosapiens"``
        - ``"musmusculus"``

    :param species:
        Species to query (see above for supported species).
        Defaults to ``"homosapiens"``.
    :type species:
        str
    :param precision:
        The level of precision to query.
        ``allele`` will query from the set of all possible alleles.
        ``gene`` will query from the set of all possible genes.
        Defaults to ``allele``.
    :type precision:
        str
    :param functionality:
        Gene/allele functionality to subset by.
        ``"any"`` queries from all possible genes/alleles.
        ``"F"`` queries from functional genes/alleles.
        ``"NF"`` queries from psuedogenes and ORFs.
        ``"P"`` queries from pseudogenes.
        ``"ORF"`` queries from ORFs.
        An allele is considered queriable if its functionality label matches the description.
        A gene is considered queriable if at least one of its alleles' functionality label matches the description.
        Defaults to ``"any"``.
    :type functionality:
        str
    :param contains_pattern:
        An optional regular expression string which will be used to filter the query result.
        If supplied, only genes/alleles which contain the regular expression will be returned.
        Defaults to ``None``.
    :type contains_pattern:
        str

    :return:
        The set of all genes/alleles that satisfy the given constraints.
    :rtype:
        FrozenSet[str]

    .. topic:: Example usage

        List all known variants for the human TR gene TRBV6-1.

        >>> tt.tr.query(species="homosapiens", contains_pattern="TRBV6-1")
        frozenset({'TRBV6-1*01'})

        List all known *Mus musculus* TRAV genes that have at least one allele which is a non-functional ORF.

        >>> tt.tr.query(species="musmusculus", precision="gene", functionality="ORF", contains_pattern="TRAV")
        frozenset({'TRAV21/DV12', 'TRAV14D-1', 'TRAV13-3', 'TRAV9D-2', 'TRAV5D-4', 'TRAV12D-3', 'TRAV12-1', 'TRAV18', 'TRAV11D'})
    """
    Parameter(species, "species").throw_error_if_not_of_type(str)
    Parameter(precision, "precision").throw_error_if_not_one_of("allele", "gene")
    Parameter(functionality, "functionality").throw_error_if_not_one_of(
        "any", "F", "NF", "P", "ORF"
    )
    Parameter(contains_pattern, "contains_pattern").throw_error_if_not_of_type(
        str, optional=True
    )

    species = _utils.clean_and_lowercase(species)

    species_is_supported = species in QUERY_ENGINES
    if not species_is_supported:
        raise ValueError(f"Unsupported species: {species}. No data available.")

    query_engine = QUERY_ENGINES[species]
    result = query_engine.query(precision, functionality)

    if contains_pattern is None:
        return result

    results_containing_substring = [i for i in result if re.search(contains_pattern, i)]
    return frozenset(results_containing_substring)
