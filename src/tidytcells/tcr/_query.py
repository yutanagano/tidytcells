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
    contains_substring: Optional[str] = None,
    contains: Optional[str] = None,
) -> FrozenSet[str]:
    """
    Query the list of all known TCR genes/alleles.

    .. topic:: Supported species

        - ``'homosapiens'``
        - ``'musmusculus'``

    :param species:
        Species to query (see above for supported species).
        Defaults to ``'homosapiens'``.
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
    :param contains_substring:
        An optional regular expression string which will be used to filter the query result.
        If supplied, only genes/alleles which contain the regular expression will be returned.
        Defaults to ``None``.
    :type contains_substring:
        str
    :param contains:
        Alias for ``contains_substring``.

        .. caution:: This will be deprecated soon in favour of ``contains_substring``.
    :type contains:
        str

    :return:
        The set of all genes/alleles that satisfy the given constraints.
    :rtype:
        FrozenSet[str]

    .. topic:: Example usage

        List all known variants for the human TCR gene TRBV6-1.

        >>> tt.tcr.query(species="homosapiens", contains_substring="TRBV6-1")
        frozenset({'TRBV6-1*01'})

        List all known *Mus musculus* TRAV genes that have at least one allele which is a non-functional ORF.

        >>> tt.tcr.query(species="musmusculus", precision="gene", functionality="ORF", contains_substring="TRAV")
        frozenset({'TRAV21/DV12', 'TRAV14D-1', 'TRAV13-3', 'TRAV9D-2', 'TRAV5D-4', 'TRAV12D-3', 'TRAV12-1', 'TRAV18', 'TRAV11D'})
    """

    contains_substring = Parameter(
        contains_substring, "contains_substring"
    ).resolve_with_alias_and_return_value(Parameter(contains, "contains"))

    Parameter(species, "species").throw_error_if_not_of_type(str)
    Parameter(precision, "precision").throw_error_if_not_one_of("allele", "gene")
    Parameter(functionality, "functionality").throw_error_if_not_one_of(
        "any", "F", "NF", "P", "ORF"
    )
    Parameter(contains_substring, "contains_substring").throw_error_if_not_of_type(
        str, optional=True
    )

    species = _utils.clean_and_lowercase(species)

    species_is_supported = species in QUERY_ENGINES
    if not species_is_supported:
        raise ValueError(f"Unsupported species: {species}. No data available.")

    query_engine = QUERY_ENGINES[species]
    result = query_engine.query(precision, functionality)

    if contains_substring is None:
        return result

    results_containing_substring = [
        i for i in result if re.search(contains_substring, i)
    ]
    return frozenset(results_containing_substring)
