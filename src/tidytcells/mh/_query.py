from functools import partial
from typing import Callable, Dict, FrozenSet, Optional, Literal

from tidytcells import config
from tidytcells import _utils
from tidytcells._utils import Parameter
from tidytcells._resources import MH_SYMBOLS_BY_SP_DB
from tidytcells.mh._query_symbols import query_symbols

_GENE_LEVEL_WARNING = (
    "tidytcells is not aware of any MH alleles for this species, and can "
    "only provide up to the level of the gene."
)
_PROTEIN_LEVEL_WARNING = (
    "tidytcells is not fully aware of all MH alleles, and the highest "
    "resolution it can provide is up to the level of the protein "
    "(two allele designations)."
)

QUERY_ENGINES: Dict[str, Callable] = {
    "homosapiens": partial(
        query_symbols, species="homosapiens", warning=_PROTEIN_LEVEL_WARNING
    ),
    "musmusculus": partial(
        query_symbols, species="musmusculus", warning=_GENE_LEVEL_WARNING
    ),
    "rattusnorvegicus": partial(
        query_symbols, species="rattusnorvegicus", warning=_GENE_LEVEL_WARNING
    ),
    "bostaurus": partial(
        query_symbols, species="bostaurus", warning=_PROTEIN_LEVEL_WARNING
    ),
    "macacamulatta": partial(
        query_symbols, species="macacamulatta", warning=_PROTEIN_LEVEL_WARNING
    ),
    "pantroglodytes": partial(
        query_symbols, species="pantroglodytes", warning=_PROTEIN_LEVEL_WARNING
    ),
    "susscrofa": partial(
        query_symbols, species="susscrofa", warning=_PROTEIN_LEVEL_WARNING
    ),
}


def query(
    species: Optional[str] = None,
    precision: Optional[Literal["allele", "gene"]] = None,
    contains_pattern: Optional[str] = None,
) -> FrozenSet[str]:
    """
    Query the list of all known MH genes / alleles.

    .. topic:: Supported species

        - ``"homosapiens"``
        - ``"bostaurus"``
        - ``"macacamulatta"``
        - ``"pantroglodytes"``
        - ``"susscrofa"``
        - ``"musmusculus"``
        - ``"rattusnorvegicus"``

    .. note::

        :py:mod:`tidytcells`' knowledge of MH alleles is limited, especially
        outside of humans. :py:mod:`tidytcells` will allow you to query HLA
        alleles up to the level of the protein (first two allele designators),
        but that is the highest resolution available. For *Mus musculus* and
        *Rattus norvegicus*, there is currently only support for gene-level
        querying.

    Parameters
    ----------
    species
        Species to query (see above for supported species). Defaults to
        ``"homosapiens"``.
    precision
        The level of precision to query. ``allele`` will query from the set of
        all possible alleles. ``gene`` will query from the set of all possible
        genes. Defaults to ``allele``.
    contains_pattern
        An optional **regular expression** string which will be used to filter
        the query result. If supplied, only genes / alleles which contain the
        regular expression will be returned. Defaults to ``None``.

    Returns
    -------
    ~typing.FrozenSet[str]
        The set of all genes / alleles that satisfy the given constraints.

    Examples
    --------

        List all known HLA-TAP1 variants (TAP is only present in the IPD
        database, which is not the default).

        >>> tt.config.set_mh_database("IPD")
        >>> tt.mh.query(species="homosapiens", contains_pattern="HLA-TAP1")
        frozenset({'HLA-TAP1*01:01', 'HLA-TAP1*01:04', 'HLA-TAP1*02:01', 'HLA-TAP1*01:06', 'HLA-TAP1*04:01', 'HLA-TAP1*03:01', 'HLA-TAP1*01:03', 'HLA-TAP1*06:01', 'HLA-TAP1*01:05', 'HLA-TAP1*01:07', 'HLA-TAP1*01:02', 'HLA-TAP1*05:01'})

        List all known *Mus musculus* haplotypes starting with H2-D.

        >>> tt.mh.query(species="musmusculus", precision="gene", contains_pattern="H2-D")
        frozenset({'H2-Dq', 'H2-Db', 'H2-Dr', 'H2-D', 'H2-Ds', 'H2-Dd', 'H2-Dk'})
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
    contains_pattern = (
        Parameter(contains_pattern, "contains_pattern")
        .throw_error_if_not_of_type(str, optional=True)
        .value
    )

    species = _utils.clean_and_lowercase(species)

    species_is_supported = species in QUERY_ENGINES
    if not species_is_supported:
        raise ValueError(f"Unsupported species: {species}. No data available.")

    mh_database = config.get_mh_database()
    if mh_database not in MH_SYMBOLS_BY_SP_DB[species]:
        raise ValueError(
            f"Unsupported mh_database: {mh_database}. No data available. "
            "Please use tidytcells.config.set_mh_database() to switch to: "
            f"{', '.join(list(MH_SYMBOLS_BY_SP_DB[species].keys()))}"
        )

    return QUERY_ENGINES[species](precision, contains_pattern, mh_database=mh_database)
