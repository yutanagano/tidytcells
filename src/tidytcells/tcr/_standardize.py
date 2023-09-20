from typing import Dict, Optional, Type

from tidytcells import _utils
from tidytcells._utils import Parameter, warnings
from tidytcells._standardized_gene_symbol import (
    StandardizedGeneSymbol,
    StandardizedHomoSapiensTrSymbol,
    StandardizedMusMusculusTrSymbol,
)


SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS: Dict[str, Type[StandardizedGeneSymbol]] = {
    "homosapiens": StandardizedHomoSapiensTrSymbol,
    "musmusculus": StandardizedMusMusculusTrSymbol,
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
        Alias for the parameter ``gene``.

        .. caution:: This will be deprecated soon in favour of ``gene``.
    :type gene_name:
        str

    :return:
        If the specified ``species`` is supported, and ``gene`` could be standardized, then return the standardized gene name.
        If ``species`` is unsupported, then the function does not attempt to standardize , and returns the unaltered ``gene`` string.
        Else follows the behaviour as set by ``on_fail``.
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

    gene = Parameter(gene, "gene").resolve_with_alias_and_return_value(
        Parameter(gene_name, "gene_name")
    )

    Parameter(gene, "gene").throw_error_if_not_of_type(str)
    Parameter(species, "species").throw_error_if_not_of_type(str)
    Parameter(enforce_functional, "enforce_functional").throw_error_if_not_of_type(bool)
    Parameter(precision, "precision").throw_error_if_not_one_of("allele", "gene")
    Parameter(on_fail, "on_fail").throw_error_if_not_one_of("reject", "keep")
    Parameter(suppress_warnings, "suppress_warnings").throw_error_if_not_of_type(bool)

    species = _utils.clean_and_lowercase(species)

    species_is_supported = species in SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS
    if not species_is_supported:
        if not suppress_warnings:
            warnings.warn_unsupported_species(species, "TCR")
        return gene

    StandardizedTrSymbolClass = SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS[species]
    standardized_tr_symbol = StandardizedTrSymbolClass(gene)

    invalid_reason = standardized_tr_symbol.get_reason_why_invalid(enforce_functional)
    if invalid_reason is not None:
        if not suppress_warnings:
            warnings.warn_failure(
                reason_for_failure=invalid_reason,
                original_input=gene,
                attempted_fix=standardized_tr_symbol.compile("allele"),
                species=species,
            )
        if on_fail == "reject":
            return None
        return gene

    return standardized_tr_symbol.compile(precision)


def standardise(*args, **kwargs):
    """
    Alias for :py:func:`tidytcells.tcr.standardize`.
    """
    return standardize(*args, **kwargs)
