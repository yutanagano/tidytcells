"""
Utility functions related to MHCs and MHC genes.
"""


import re
from .._resources import *
from typing import Optional, FrozenSet
from .._utils.abstract_functions import standardize_template, query_template
from .._utils.gene_query_engines import HLAQueryEngine, MusMusculusMHCQueryEngine
from .._utils.gene_standardizers import HLAStandardizer, MusMusculusMHCStandardizer
from .._utils.warnings import *


# --- STATIC RESOURCES ---

CHAIN_ALPHA_RE = re.compile(r"HLA-([ABCEFG]|D[PQR]A)")
CHAIN_BETA_RE = re.compile(r"HLA-D[PQR]B|B2M")
CLASS_1_RE = re.compile(r"HLA-[ABCEFG]|B2M")
CLASS_2_RE = re.compile(r"HLA-D[PQR][AB]")


# --- MAIN FUNCTIONS ---

STANDARDIZERS = {
    "homosapiens": HLAStandardizer,
    "musmusculus": MusMusculusMHCStandardizer,
}

QUERY_ENGINES = {
    "homosapiens": HLAQueryEngine,
    "musmusculus": MusMusculusMHCQueryEngine,
}


def standardize(
    gene: Optional[str] = None,
    species: str = "homosapiens",
    precision: str = "allele",
    on_fail: str = "reject",
    suppress_warnings: bool = False,
    gene_name: Optional[str] = None,
) -> tuple:
    """
    Attempt to standardize an MHC gene name to be IMGT-compliant.

    .. topic:: Supported species

        - ``'homosapiens'``
        - ``'musmusculus'``

    .. note::
        This function will only verify the validity of an MHC gene/allele up to the level of the protein.
        Any further precise allele designations will not be verified, apart from the requirement that the format (colon-separated numbers) look valid.
        The reasons for this is firstly because new alleles at that level are added to the IMGT list quite often and so accurate verification is difficult, secondly because people rarely need verification to such a precise level, and finally because such verification costs more computational effort with diminishing returns.

    :param gene:
        Potentially non-standardized MHC gene name.
    :type gene:
        str
    :param species:
        Species to which the MHC gene belongs (see above for supported species).
        Defaults to ``'homosapiens'``.
    :type species:
        str
    :param precision:
        The maximum level of precision to standardize to.
        ``'allele'`` standardizes to the maximum precision possible.
        ``'protein'`` keeps allele designators up to the level of the protein.
        ``'gene'`` standardizes only to the level of the gene.
        Defaults to ``'allele'``.
    :type precision:
        str
    :param on_fail:
        Behaviour when standardization fails.
        If set to ``"reject"``, returns ``None`` on failure.
        If set to ``"keep"``, returns the original input.
        Defaults to ``"reject"``.
    :type on_fail:
        str
    :param suppress_warnings:
        Disable warnings that are usually emitted when standardisation fails.
        Defaults to ``False``.
    :type suppress_warnings:
        bool

    :param gene_name:
        Alias for the parameter ``gene``.

        .. caution:: This will be deprecated soon in favour of ``gene``.
    :type gene_name:
        str

    :return:
        If the specified ``species`` is supported, and ``gene`` could be standardized, then return the standardized gene name.
        If ``species`` is unsupported, then the function does not attempt to standardize, and returns the unaltered ``gene`` string.
        Else follows the behvaiour as set by ``on_fail``.
    :rtype:
        Union[str, None]

    .. topic:: Example usage

        Input strings will intelligently be corrected to IMGT-compliant gene symbols.

        >>> tt.mhc.standardize("A1")
        'HLA-A*01'

        The ``precision`` setting can truncate unnecessary information.

        >>> tt.mhc.standardize("HLA-A*01", precision="gene")
        'HLA-A'

        *Mus musculus* is a supported species.

        >>> tt.mhc.standardize("CRW2", species="musmusculus")
        'MH1-M5'
    """
    # Alias resolution
    if gene is None and not gene_name is None:
        warn(
            'The parameter "gene_name" will be deprecated in the near future. Please switch to using "gene".',
            FutureWarning,
        )
        gene = gene_name

    return standardize_template(
        gene=gene,
        gene_type="MHC",
        species=species,
        enforce_functional=True,
        precision=precision,
        on_fail=on_fail,
        suppress_warnings=suppress_warnings,
        standardizer_dict=STANDARDIZERS,
        allowed_precision={"allele", "protein", "gene"},
    )


def query(
    species: str = "homosapiens",
    precision: str = "allele",
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
    :param contains:
        An optional **regular expression** string which will be used to filter the query result.
        If supplied, only genes/alleles which contain the regular expression will be returned.
        Defaults to ``None``.
    :type contains:
        str

    :return:
        The set of all genes/alleles that satisfy the given constraints.
    :rtype:
        FrozenSet[str]

    .. topic:: Example usage

        List all known HLA-G variants.

        >>> tt.mhc.query(species="homosapiens", contains="HLA-G")
        frozenset({'HLA-TAP1*03:01', 'HLA-TAP1*01:02', 'HLA-TAP1*06:01', 'HLA-TAP1*04:01', 'HLA-TAP1*02:01', 'HLA-TAP1*05:01', 'HLA-TAP1*01:01'})

        List all known *Mus musculus* MH1-Q genes.

        >>> tt.mhc.query(species="musmusculus", precision="gene", contains="MH1-Q")
        frozenset({'MH1-Q3', 'MH1-Q9', 'MH1-Q1', 'MH1-Q2', 'MH1-Q6', 'MH1-Q10', 'MH1-Q5', 'MH1-Q8', 'MH1-Q7', 'MH1-Q4'})
    """

    return query_template(
        species=species,
        precision=precision,
        functionality="any",
        contains=contains,
        query_engine_dict=QUERY_ENGINES,
    )


def get_chain(
    gene: Optional[str] = None,
    suppress_warnings: bool = False,
    gene_name: Optional[str] = None,
) -> str:
    """
    Given a standardized MHC gene name, detect whether it codes for an alpha or a beta chain molecule.

    .. note::

        This function currently only recognises HLAs, and not MHCs from other species.

    :param gene:
        Standardized MHC gene name
    :type gene:
        str
    :param suppress_warnings:
        Disable warnings that are usually emitted when chain classification fails.
        Defaults to ``False``.
    :type suppress_warnings:
        bool

    :param gene_name:
        Alias for the parameter ``gene``.

        .. caution:: This will be deprecated soon in favour of ``gene``.
    :type gene_name:
        str

    :return:
        ``'alpha'`` or ``'beta'`` if ``gene`` is recognised and its chain is known, else ``None``.
    :rtype:
        Union[str, None]

    .. topic:: Example usage

        >>> tt.mhc.get_chain("HLA-A")
        'alpha'
        >>> tt.mhc.get_chain("HLA-DRB2")
        'beta'
        >>> tt.mhc.get_chain("B2M")
        'beta'
    """
    # Alias resolution
    if gene is None and not gene_name is None:
        warn(
            'The parameter "gene_name" will be deprecated in the near future. Please switch to using "gene".',
            FutureWarning,
        )
        gene = gene_name

    if type(gene) == str:
        # If we don't recognise the gene, return None with warning
        gene = gene.split("*")[0]

        if not gene in (*HOMOSAPIENS_MHC, "B2M"):
            if not suppress_warnings:
                warn(f"Unrecognised gene {gene}. Is this standardized?")
            return None

        if CHAIN_ALPHA_RE.match(gene):
            return "alpha"

        if CHAIN_BETA_RE.match(gene):
            return "beta"

        if not suppress_warnings:
            warn(f"Chain for {gene} unknown.")
        return None

    raise TypeError(f"gene_name must be type str, got {gene} ({type(gene)}).")


def get_class(
    gene: Optional[str] = None,
    suppress_warnings: bool = False,
    gene_name: Optional[str] = None,
) -> int:
    """
    Given a standardized MHC gene name, detect whether it comprises a class I
    or II MHC receptor complex.

    .. note::

        This function currently only recognises HLAs, and not MHCs from other species.

    :param gene:
        Standardized MHC gene name
    :type gene:
        str
    :param suppress_warnings:
        Disable warnings that are usually emitted when classification fails.
        Defaults to ``False``.
    :type suppress_warnings:
        bool

    :param gene_name:
        Alias for the parameter ``gene``.

        .. caution:: This will be deprecated soon in favour of ``gene``.
    :type gene_name:
        str

    :return:
        ``1`` or ``2`` if ``gene`` is recognised and its class is known, else ``None``.
    :rtype:
        Union[int, None]

    .. topic:: Example usage

        >>> tt.mhc.get_class("HLA-A")
        1
        >>> tt.mhc.get_class("HLA-DRB2")
        2
        >>> tt.mhc.get_class("B2M")
        1
    """
    # Alias resolution
    if gene is None and not gene_name is None:
        warn(
            'The parameter "gene_name" will be deprecated in the near future. Please switch to using "gene".',
            FutureWarning,
        )
        gene = gene_name

    if type(gene) == str:
        # If we don't recognise the gene, return None with warning
        gene = gene.split("*")[0]

        if not gene in (*HOMOSAPIENS_MHC, "B2M"):
            if not suppress_warnings:
                warn(f"Unrecognised gene {gene}. Is this standardized?")
            return None

        if CLASS_1_RE.match(gene):
            return 1

        if CLASS_2_RE.match(gene):
            return 2

        if not suppress_warnings:
            warn(f"Class for {gene} unknown.")
        return None

    raise TypeError(f"gene_name must be type str, got {gene} ({type(gene)}).")
