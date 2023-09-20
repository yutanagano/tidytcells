import re
from typing import Dict, FrozenSet, Optional, Type

from tidytcells import _utils
from tidytcells._utils import Parameter
from tidytcells._utils.gene_query_engines import (
    GeneQueryEngine,
    HLAQueryEngine,
    MusMusculusMHCQueryEngine,
)


QUERY_ENGINES: Dict[str, Type[GeneQueryEngine]] = {
    "homosapiens": HLAQueryEngine,
    "musmusculus": MusMusculusMHCQueryEngine,
}


def query(
    species: str = "homosapiens",
    precision: str = "allele",
    contains_substring: Optional[str] = None,
    contains: Optional[str] = None,
) -> FrozenSet[str]:
    """
    Query the list of all known MHC genes/alleles.

    .. topic:: Supported species

        - ``'homosapiens'``
        - ``'musmusculus'``

    .. note::

        :py:mod:`tidytcells`' knowledge of MHC alleles is limited, especially outside of humans.
        :py:mod:`tidytcells` will allow you to query HLA alleles up to the level of the protein (first two allele designators), but that is the highest resolution available.
        For Mus musculus, there is currently only support for gene-level querying.

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
    :param contains_substring:
        An optional **regular expression** string which will be used to filter the query result.
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

        List all known HLA-G variants.

        >>> tt.mhc.query(species="homosapiens", contains_substring="HLA-G")
        frozenset({'HLA-TAP1*03:01', 'HLA-TAP1*01:02', 'HLA-TAP1*06:01', 'HLA-TAP1*04:01', 'HLA-TAP1*02:01', 'HLA-TAP1*05:01', 'HLA-TAP1*01:01'})

        List all known *Mus musculus* MH1-Q genes.

        >>> tt.mhc.query(species="musmusculus", precision="gene", contains_substring="MH1-Q")
        frozenset({'MH1-Q3', 'MH1-Q9', 'MH1-Q1', 'MH1-Q2', 'MH1-Q6', 'MH1-Q10', 'MH1-Q5', 'MH1-Q8', 'MH1-Q7', 'MH1-Q4'})
    """
    contains_substring = Parameter(
        contains_substring, "contains_substring"
    ).resolve_with_alias_and_return_value(Parameter(contains, "contains"))

    Parameter(species, "species").throw_error_if_not_of_type(str)
    Parameter(precision, "precision").throw_error_if_not_one_of("allele", "gene")
    Parameter(contains_substring, "contains_substring").throw_error_if_not_of_type(
        str, optional=True
    )

    species = _utils.lowercase_and_remove_whitespace(species)

    species_is_supported = species in QUERY_ENGINES
    if not species_is_supported:
        raise ValueError(f"Unsupported species: {species}. No data available.")

    query_engine = QUERY_ENGINES[species]
    result = query_engine.query(precision, functionality=None)

    if contains_substring is None:
        return result

    results_containing_substring = [
        i for i in result if re.search(contains_substring, i)
    ]
    return frozenset(results_containing_substring)
