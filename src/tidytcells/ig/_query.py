import re
from typing import Dict, FrozenSet, Optional, Type

from tidytcells import _utils
from tidytcells._utils import Parameter
from tidytcells._query_engine import (
    QueryEngine,
    HomoSapiensIgQueryEngine,
)


QUERY_ENGINES: Dict[str, Type[QueryEngine]] = {
    "homosapiens": HomoSapiensIgQueryEngine,
}


def query(
    species: Optional[str] = None,
    precision: Optional[str] = None,
    functionality: Optional[str] = None,
    contains_pattern: Optional[str] = None,
) -> FrozenSet[str]:
    """
    Query the list of all known IG
     genes / alleles.

    .. topic:: Supported species

        - ``"homosapiens"``

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
        The set of all genes / alleles that satisfy the given constraints.
    :rtype:
        FrozenSet[str]

    .. topic:: Example usage

        List all known variants for the human IG gene IGLJ1.

        >>> tt.ig.query(species="homosapiens", contains_pattern="IGLJ1")
        frozenset({'IGL1*01'})

        List all known IGLV genes that have at least one allele which is a non-functional ORF.

        >>> tt.ig.query(species="homosapiens", precision="gene", functionality="ORF", contains_pattern="IGLV")
        frozenset({'IGLV1-41', 'IGLV1-50', 'IGLV11-55', 'IGLV2-33', 'IGLV3-32', 'IGLV5-48', 'IGLV8/OR8-1'})
    """
    species = (
        Parameter(species, "species")
        .set_default("homosapiens")
        .throw_error_if_not_of_type(str)
        .value
    )
    precision = (
        Parameter(precision, "precision")
        .set_default("allele")
        .throw_error_if_not_one_of("allele", "gene")
        .value
    )
    functionality = (
        Parameter(functionality, "functionality")
        .set_default("any")
        .throw_error_if_not_one_of("any", "F", "NF", "P", "ORF")
        .value
    )
    contains_pattern = (
        Parameter(contains_pattern, "contains_pattern")
        .throw_error_if_not_of_type(str, optional=True)
        .value
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
