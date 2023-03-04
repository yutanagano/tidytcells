"""
Utility functions related to TCRs and TCR genes.
"""


from typing import FrozenSet, Optional
from .._utils.abstract_functions import standardise_template, query_template
from .._utils.gene_query_engines import (
    HomoSapiensTCRQueryEngine,
    MusMusculusTCRQueryEngine,
)
from .._utils.gene_standardisers import (
    HomoSapiensTCRStandardiser,
    MusMusculusTCRStandardiser,
)
from .._utils.warnings import *


STANDARDISERS = {
    "homosapiens": HomoSapiensTCRStandardiser,
    "musmusculus": MusMusculusTCRStandardiser,
}

QUERY_ENGINES = {
    "homosapiens": HomoSapiensTCRQueryEngine,
    "musmusculus": MusMusculusTCRQueryEngine,
}


def standardise(
    gene: Optional[str] = None,
    species: str = "homosapiens",
    enforce_functional: bool = False,
    precision: str = "allele",
    suppress_warnings: bool = False,
    gene_name: Optional[str] = None,
) -> str:
    """
    Attempt to standardise a TCR gene name to be IMGT-compliant.

    .. topic:: Supported species

        - ``'homosapiens'``
        - ``'musmusculus'``

    :param gene:
        Potentially non-standardised TCR gene name.
    :type gene:
        ``str``
    :param species:
        Species to which the TCR gene belongs (see above for supported species).
        Defaults to ``'homosapiens'``.
    :type species:
        ``str``
    :param enforce_functional:
        If ``True``, disallows TCR genes that are recognised by IMGT but are marked as non-functional (ORF or pseudogene).
        Defaults to ``False``.
    :type enforce_functional:
        ``bool``
    :param precision:
        The maximum level of precision to standardise to.
        ``'allele'`` standardises to the maximum precision possible.
        ``'gene'`` standardises only to the level of the gene.
        Defaults to ``'allele'``.
    :type precision:
        ``str``
    :param suppress_warnings:
        Disable warnings that are usually emitted when standardisation fails.
        Defaults to ``False``.
    :type suppress_warnings:
        ``bool``

    :param gene_name:
        Alias for the parameter ``gene``.
    :type gene_name:
        ``str``

    :return:
        If the specified ``species`` is supported, and ``gene`` could be standardised, then return the standardised gene name.
        If ``species`` is unsupported, then the function does not attempt to standardise , and returns the unaltered ``gene`` string.
        Else returns ``None``.
    :rtype:
        ``str`` or ``None``
    """
    # Alias resolution
    if gene is None:
        gene = gene_name

    return standardise_template(
        gene=gene,
        gene_type="TCR",
        species=species,
        enforce_functional=enforce_functional,
        precision=precision,
        suppress_warnings=suppress_warnings,
        standardiser_dict=STANDARDISERS,
        allowed_precision={"allele", "gene"},
    )


def query(
    species: str = "homosapiens",
    precision: str = "allele",
    functionality: str = "any",
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
        ``str``
    :param precision:
        The level of precision to query.
        ``allele`` will query from the set of all possible alleles.
        ``gene`` will query from the set of all possible genes.
        Defaults to ``allele``.
    :type precision:
        ``str``
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
        ``str``
    :param contains:
        An optional regular expression string which will be used to filter the query result.
        If supplied, only genes/alleles which contain the regular expression will be returned.
        Defaults to ``None``.
    :type contains:
        ``str``

    :return:
        The set of all genes/alleles that satisfy the given constraints.
    :rtype:
        ``FrozenSet[str]``
    """

    return query_template(
        species=species,
        precision=precision,
        functionality=functionality,
        contains=contains,
        query_engine_dict=QUERY_ENGINES,
    )
