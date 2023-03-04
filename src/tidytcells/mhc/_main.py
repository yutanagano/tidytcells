"""
Utility functions related to MHCs and MHC genes.
"""


import re
from .._resources import *
from typing import Optional, FrozenSet
from .._utils.abstract_functions import standardise_template, query_template
from .._utils.gene_query_engines import HLAQueryEngine, MusMusculusMHCQueryEngine
from .._utils.gene_standardisers import HLAStandardiser, MusMusculusMHCStandardiser
from .._utils.warnings import *


# --- STATIC RESOURCES ---

CHAIN_ALPHA_RE = re.compile(r"HLA-([ABCEFG]|D[PQR]A)")
CHAIN_BETA_RE = re.compile(r"HLA-D[PQR]B|B2M")
CLASS_1_RE = re.compile(r"HLA-[ABCEFG]|B2M")
CLASS_2_RE = re.compile(r"HLA-D[PQR][AB]")


# --- MAIN FUNCTIONS ---

STANDARDISERS = {
    "homosapiens": HLAStandardiser,
    "musmusculus": MusMusculusMHCStandardiser,
}

QUERY_ENGINES = {
    "homosapiens": HLAQueryEngine,
    "musmusculus": MusMusculusMHCQueryEngine
}


def standardise(
    gene: Optional[str] = None,
    species: str = "homosapiens",
    precision: str = "allele",
    suppress_warnings: bool = False,
    gene_name: Optional[str] = None,
) -> tuple:
    """
    Attempt to standardise an MHC gene name to be IMGT-compliant.

    .. note::
        This function will only verify the validity of an MHC gene/allele up to the level of the protein.
        Any further precise allele designations will not be verified, apart from the requirement that the format (colon-separated numbers) look valid.
        The reasons for this is firstly because new alleles at that level are added to the IMGT list quite often and so accurate verification is difficult, secondly because people rarely need verification to such a precise level, and finally because such verification costs more computational effort with diminishing returns.

    :param gene:
        Potentially non-standardised MHC gene name.
    :type gene:
        ``str``
    :param species:
        Species to which the MHC gene belongs (see :ref:`supported_species`).
        Defaults to ``'homosapiens'``.
    :type species:
        ``str``
    :param precision:
        The maximum level of precision to standardise to.
        ``'allele'`` standardises to the maximum precision possible.
        ``'protein'`` keeps allele designators up to the level of the protein.
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
        If ``species`` is unsupported, then the function does not attempt to standardise, and returns the unaltered ``gene`` string.
        Else returns ``None``.
    :rtype:
        ``str`` or ``None``
    """
    # Alias resolution
    if gene is None:
        gene = gene_name

    return standardise_template(
        gene=gene,
        gene_type="MHC",
        species=species,
        enforce_functional=True,
        precision=precision,
        suppress_warnings=suppress_warnings,
        standardiser_dict=STANDARDISERS,
        allowed_precision={'allele', 'protein', 'gene'}
    )


def query(species: str = "homosapiens", precision: str = "protein") -> FrozenSet[str]:
    """
    Query the list of all known MHC genes/alleles.

    :param species:
        Species to query (see :ref:`supported_species`).
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

    :return:
        The set of all genes/alleles that satisfy the given constraints.
    :rtype:
        ``FrozenSet[str]``
    """

    return query_template(species=species, precision=precision, query_engine_dict=QUERY_ENGINES)


def get_chain(
    gene: Optional[str] = None,
    suppress_warnings: bool = False,
    gene_name: Optional[str] = None,
) -> str:
    """
    Given a standardised MHC gene name, detect whether it codes for an alpha or a beta chain molecule.

    :param gene:
        Standardised MHC gene name
    :type gene:
        ``str``
    :param suppress_warnings:
        Disable warnings that are usually emitted when chain classification fails.
        Defaults to ``False``.
    :type suppress_warnings:
        ``bool``

    :param gene_name:
        Alias for the parameter ``gene``.
    :type gene_name:
        ``str``

    :return:
        ``'alpha'`` or ``'beta'`` if ``gene`` is recognised and its chain is known, else ``None``.
    :rtype:
        ``str`` or ``None``
    """
    # Alias resolution
    if gene is None:
        gene = gene_name

    if type(gene) == str:
        # If we don't recognise the gene, return None with warning
        gene = gene.split("*")[0]

        if not gene in (*HOMOSAPIENS_MHC, "B2M"):
            if not suppress_warnings:
                warn(f"Unrecognised gene {gene}. Is this standardised?")
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
    Given a standardised MHC gene name, detect whether it comprises a class I
    or II MHC receptor complex.

    :param gene:
        Standardised MHC gene name
    :type gene:
        ``str``
    :param suppress_warnings:
        Disable warnings that are usually emitted when classification fails.
        Defaults to ``False``.
    :type suppress_warnings:
        ``bool``

    :param gene_name:
        Alias for the parameter ``gene``.
    :type gene_name:
        ``str``

    :return:
        ``1`` or ``2`` if ``gene`` is recognised and its class is known, else ``None``.
    :rtype:
        ``int`` or ``None``
    """
    # Alias resolution
    if gene is None:
        gene = gene_name

    if type(gene) == str:
        # If we don't recognise the gene, return None with warning
        gene = gene.split("*")[0]

        if not gene in (*HOMOSAPIENS_MHC, "B2M"):
            if not suppress_warnings:
                warn(f"Unrecognised gene {gene}. Is this standardised?")
            return None

        if CLASS_1_RE.match(gene):
            return 1

        if CLASS_2_RE.match(gene):
            return 2

        if not suppress_warnings:
            warn(f"Class for {gene} unknown.")
        return None

    raise TypeError(f"gene_name must be type str, got {gene} ({type(gene)}).")
