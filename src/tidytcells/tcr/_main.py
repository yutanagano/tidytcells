"""
Utility functions related to TCRs and TCR genes.
"""


from typing import Dict, FrozenSet, Optional
from .._resources import HOMOSAPIENS_TCR_AA_SEQUENCES
from .._utils.abstract_functions import standardize_template, query_template
from .._utils.gene_query_engines import (
    HomoSapiensTCRQueryEngine,
    MusMusculusTCRQueryEngine,
)
from .._utils.gene_standardizers import (
    HomoSapiensTCRStandardizer,
    MusMusculusTCRStandardizer,
)
from .._utils.warnings import *


STANDARDIZERS = {
    "homosapiens": HomoSapiensTCRStandardizer,
    "musmusculus": MusMusculusTCRStandardizer,
}

QUERY_ENGINES = {
    "homosapiens": HomoSapiensTCRQueryEngine,
    "musmusculus": MusMusculusTCRQueryEngine,
}


def standardize(
    gene: Optional[str] = None,
    species: str = "homosapiens",
    enforce_functional: bool = False,
    precision: str = "allele",
    on_fail: str = "reject",
    suppress_warnings: bool = False,
    gene_name: Optional[str] = None,
) -> str:
    """
    Attempt to standardize a TCR gene name to be IMGT-compliant.

    .. topic:: Supported species

        - ``'homosapiens'``
        - ``'musmusculus'``

    :param gene:
        Potentially non-standardized TCR gene name.
    :type gene:
        str
    :param species:
        Species to which the TCR gene belongs (see above for supported species).
        Defaults to ``'homosapiens'``.
    :type species:
        str
    :param enforce_functional:
        If ``True``, disallows TCR genes that are recognised by IMGT but are marked as non-functional (ORF or pseudogene).
        Defaults to ``False``.
    :type enforce_functional:
        bool
    :param precision:
        The maximum level of precision to standardize to.
        ``'allele'`` standardizes to the maximum precision possible.
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
        Alias for the parameter ``gene``. This will be deprecated soon.
    :type gene_name:
        str

    :return:
        If the specified ``species`` is supported, and ``gene`` could be standardized, then return the standardized gene name.
        If ``species`` is unsupported, then the function does not attempt to standardize , and returns the unaltered ``gene`` string.
        Else returns ``None``.
    :rtype:
        Union[str, None]

    .. topic:: Example usage

        Input strings will intelligently be corrected to IMGT-compliant gene symbols.

        >>> tt.tcr.standardize("aj1")
        'TRAJ1'

        The ``precision`` setting can truncate unnecessary information.

        >>> tt.tcr.standardize("TRBV6-4*01", precision="gene")
        'TRBV6-4'

        The ``enforce_functional`` setting will cause non-functional genes or alleles to be rejected.

        >>> result = tt.tcr.standardize("TRBV1", enforce_functional=True)
        UserWarning: Failed to standardize "TRBV1" for species homosapiens: gene has no functional alleles. Attempted fix "TRBV1".
        >>> print(result)
        None

        *Mus musculus* is a supported species.

        >>> tt.tcr.standardize("TCRBV22S1A2N1T", species="musmusculus")
        'TRBV2'
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
        gene_type="TCR",
        species=species,
        enforce_functional=enforce_functional,
        precision=precision,
        on_fail=on_fail,
        suppress_warnings=suppress_warnings,
        standardizer_dict=STANDARDIZERS,
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
    :param contains:
        An optional regular expression string which will be used to filter the query result.
        If supplied, only genes/alleles which contain the regular expression will be returned.
        Defaults to ``None``.
    :type contains:
        str

    :return:
        The set of all genes/alleles that satisfy the given constraints.
    :rtype:
        FrozenSet[str]

    .. topic:: Example usage

        List all known variants for the human TCR gene TRBV6-1.

        >>> tt.tcr.query(species="homosapiens", contains="TRBV6-1")
        frozenset({'TRBV6-1*01'})

        List all known *Mus musculus* TRAV genes that have at least one allele which is a non-functional ORF.

        >>> tt.tcr.query(species="musmusculus", precision="gene", functionality="ORF", contains="TRAV")
        frozenset({'TRAV21/DV12', 'TRAV14D-1', 'TRAV13-3', 'TRAV9D-2', 'TRAV5D-4', 'TRAV12D-3', 'TRAV12-1', 'TRAV18', 'TRAV11D'})
    """

    return query_template(
        species=species,
        precision=precision,
        functionality=functionality,
        contains=contains,
        query_engine_dict=QUERY_ENGINES,
    )


def get_aa_sequence(gene: str, species: str = "homosapiens") -> Dict[str, str]:
    """
    Look up the amino acid sequence of a given TCR gene.

    .. topic:: Supported species

        - ``'homosapiens'``

    .. note::

        This function currently only supports V genes.
        Support for J genes is planned for the future.

    :param gene:
        Standardized gene name.
        The gene must be specified to the level of the allele.
        Note that some genes, notably the non-functional ones, will not have resolvable amino acid sequences.
    :type gene:
        str
    :param species:
        Species to which the TCR gene in question belongs (see above for supported species).
        Defaults to ``'homosapiens'``.
    :type species:
        str

    :return:
        A dictionary with keys corresponding to names of different sequence regions within the gene, and values corresponding to their amino acid sequences.
    :rtype:
        Dict[str, str]

    .. topic:: Example usage

        Get amino acid sequence information about the human V gene TRBV2*01.

        >>> tt.tcr.get_aa_sequence(gene="TRBV2*01", species="homosapiens")
        {'CDR1-IMGT': 'SNHLY', 'CDR2-IMGT': 'FYNNEI', 'FR1-IMGT': 'EPEVTQTPSHQVTQMGQEVILRCVPI', 'FR2-IMGT': 'FYWYRQILGQKVEFLVS', 'FR3-IMGT': 'SEKSEIFDDQFSVERPDGSNFTLKIRSTKLEDSAMYFC', 'V-REGION': 'EPEVTQTPSHQVTQMGQEVILRCVPISNHLYFYWYRQILGQKVEFLVSFYNNEISEKSEIFDDQFSVERPDGSNFTLKIRSTKLEDSAMYFCASSE'}
    """
    # Type checks
    if type(gene) != str:
        raise TypeError(f"gene must be type str, got {gene} ({type(gene)}).")
    if type(species) != str:
        raise TypeError(f"species must be type str, got {species} ({type(species)}).")

    # For backward compatibility, fix CamelCased species
    species = "".join(species.split()).lower()

    # Currently only supports homosapiens
    if species != "homosapiens":
        raise ValueError(f"Unsupported species: {species}. No data available.")

    if gene in HOMOSAPIENS_TCR_AA_SEQUENCES:
        return HOMOSAPIENS_TCR_AA_SEQUENCES[gene]

    raise ValueError(f"No data found for TCR gene {gene} for species {species}.")
