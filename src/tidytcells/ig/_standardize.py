import logging
from tidytcells import _utils
from tidytcells._utils import Parameter
from tidytcells._standardized_gene_symbol import (
    StandardizedSymbol,
    StandardizedHomoSapiensIgSymbol,
)
from typing import Dict, Optional, Type, Literal


logger = logging.getLogger(__name__)


SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS: Dict[str, Type[StandardizedSymbol]] = {
    "homosapiens": StandardizedHomoSapiensIgSymbol,
}


def standardize(
    symbol: Optional[str] = None,
    species: Optional[str] = None,
    enforce_functional: Optional[bool] = None,
    allow_subgroup: Optional[bool] = None,
    precision: Optional[Literal["allele", "gene", "subgroup"]] = None,
    on_fail: Optional[Literal["reject", "keep"]] = None,
    log_failures: Optional[bool] = None,
    gene: Optional[str] = None,
    suppress_warnings: Optional[bool] = None,
) -> Optional[str]:
    """
    Attempt to standardize a IG gene / allele symbol to be IMGT-compliant.

    .. topic:: Supported species

        - ``"homosapiens"``

    :param symbol:
        Potentially non-standardized IG gene / allele symbol.
    :type symbol:
        str
    :param species:
        Can be specified to standardise to a IG symbol that is known to be valid for that species (see above for supported species).
        Currently, only *Homo sapiens* is supported, but this parameter has been kept to keep the interface compatible with that of its sister function in :py:mod:`tidytcells.tr`.
        Defaults to ``"homosapiens"``.
    :type species:
        str
    :param enforce_functional:
        If ``True``, disallows IG genes / alleles that are recognised by IMGT but are marked as non-functional (ORF or pseudogene).
        Defaults to ``False``.
    :type enforce_functional:
        bool
    :param allow_subgroup:
        If ``True``, allows valid subgroups (as well as more specific gene/allele symbos) to pass standardisation.
        If ``False``, the supplied symbol must point to at least a specific gene.
        Defaults to ``False``.
    :type allow_subgroup:
        bool
    :param precision:
        The maximum level of precision to standardize to.
        ``"allele"`` standardizes to the maximum precision possible.
        ``"gene"`` standardizes only to the level of the gene.
        Defaults to ``"allele"``.
    :type precision:
        str
    :param on_fail:
        Behaviour when standardization fails.
        If set to ``"reject"``, returns ``None`` on failure.
        If set to ``"keep"``, returns the original input.
        Defaults to ``"reject"``.
    :type on_fail:
        str
    :param log_failures:
        Report standardisation failures through logging (at level ``WARNING``).
        Defaults to ``True``.
    :type log_failures:
        bool
    :param gene:
        Alias for the parameter `symbol`.
    :type gene:
        str
    :param suppress_warnings:
        Disable warnings that are usually logged when standardisation fails.
        Deprecated in favour of `log_failures`.
    :type suppress_warnings:
        bool

    :return:
        If the specified `species` is supported, and `symbol` could be standardized, then return the standardized symbol name.
        If `species` is unsupported, then the function does not attempt to standardize , and returns the unaltered `symbol` string.
        Else follows the behaviour as set by `on_fail`.
    :rtype:
        Optional[str]

    .. topic:: Example usage

        Input strings will intelligently be corrected to IMGT-compliant gene / allele symbols.

        >>> tt.ig.standardize("lj1")
        'IGLJ1'

        The `precision` setting can truncate unnecessary information.

        >>> tt.ig.standardize("IGHV1-18*02", precision="gene")
        'IGHV1-18'

        The `enforce_functional` setting will cause non-functional genes or alleles to be rejected.

        >>> result = tt.ig.standardize("IGHV1-12", enforce_functional=True)
        Failed to standardize "IGHV1-12" for species homosapiens: gene has no functional alleles. Attempted fix "IGHV1-12".
        >>> print(result)
        None

    .. topic:: Decision Logic

        To provide an easy way to gauge the scope and limitations of standardization, below is a simplified overview of the decision logic employed when attempting to standardize a TR symbol.
        For more detail, please refer to the `source code <https://github.com/yutanagano/tidytcells>`_.

        .. code-block:: none

            IF the specified species is not supported for standardization:
                RETURN original symbol without modification

            ELSE:
                // attempt standardization
                {
                    IF symbol is already in IMGT-compliant form:
                        set standardization status as successful
                        skip rest of standardization

                    IF symbol is a known deprecated symbol:
                        overwrite symbol with current IMGT-compliant symbol
                        set standardization status as successful
                        skip rest of standardization

                    replace "." with "-"                                        //e.g. IGHV1.2 -> IGHV1-2
                    add back any missing backslashes                            //e.g. IGHV1OR15-1 -> IGHV1/OR15-1
                    remove any unnecessary trailing zeros                       //e.g. IGHV1-02 -> IGHV1-2
                    IF symbol is now in IMGT-compliant form:
                        set standardization status as successful
                        skip rest of standardization


                    add "IG" to the beginning of the symbol if necessary   //e.g. HV1-18 -> IGHV1-18
                    IF symbol is now in IMGT-compliant form:
                        set standardization status as successful
                        skip rest of standardization

                    try adding or removing "-1" from the end of the symbol //e.g. IGHV6 -> IGHV6-1
                    IF symbol is now in IMGT-compliant form:
                        set standardization status as successful
                        skip rest of standardization

                    set standardization status as failed
                }

                IF standardization status is set to successful:
                    RETURN standardized symbol

                ELSE:
                    IF on_fail is set to "reject":
                        RETURN None
                    IF on_fail is set to "keep":
                        RETURN original symbol without modification
    """
    symbol = (
        Parameter(symbol, "symbol")
        .resolve_with_alias(gene, "gene")
        .throw_error_if_not_of_type(str)
        .value
    )
    species = (
        Parameter(species, "species")
        .set_default("homosapiens")
        .throw_error_if_not_of_type(str)
        .value
    )
    enforce_functional = (
        Parameter(enforce_functional, "enforce_functional")
        .set_default(False)
        .throw_error_if_not_of_type(bool)
        .value
    )
    allow_subgroup = (
        Parameter(allow_subgroup, "allow_subgroup")
        .set_default(False)
        .throw_error_if_not_of_type(bool)
        .value
    )
    precision = (
        Parameter(precision, "precision")
        .set_default("allele")
        .throw_error_if_not_one_of("allele", "gene", "subgroup")
        .value
    )
    on_fail = (
        Parameter(on_fail, "on_fail")
        .set_default("reject")
        .throw_error_if_not_one_of("reject", "keep")
        .value
    )
    suppress_warnings_inverted = (
        not suppress_warnings if suppress_warnings is not None else None
    )
    log_failures = (
        Parameter(log_failures, "log_failures")
        .set_default(True)
        .resolve_with_alias(suppress_warnings_inverted, "suppress_warnings")
        .throw_error_if_not_of_type(bool)
        .value
    )

    allow_subgroup = True if precision == "subgroup" else allow_subgroup
    species = _utils.clean_and_lowercase(species)

    if species == "any":
        best_attempt_invalid_reason = None
        best_attempt_standardised_symbol = None
        best_attempt_species = None

        for (
            species,
            StandardizedIgSymbolClass,
        ) in SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS.items():
            standardized_ig_symbol = StandardizedIgSymbolClass(symbol, allow_subgroup)
            invalid_reason = standardized_ig_symbol.get_reason_why_invalid(
                enforce_functional
            )

            if invalid_reason is None:
                return standardized_ig_symbol.compile(precision)

            if species == "homosapiens":
                best_attempt_invalid_reason = invalid_reason
                best_attempt_standardised_symbol = standardized_ig_symbol
                best_attempt_species = species

        if log_failures:
            _utils.warn_failure(
                reason_for_failure=best_attempt_invalid_reason,
                original_input=symbol,
                attempted_fix=best_attempt_standardised_symbol.compile("allele"),
                species=best_attempt_species,
                logger=logger,
            )
        if on_fail == "reject":
            return None
        return symbol

    if species not in SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS:
        if log_failures:
            _utils.warn_unsupported_species(species, "IG", logger)
        return symbol

    StandardizedIgSymbolClass = SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS[species]
    standardized_ig_symbol = StandardizedIgSymbolClass(symbol, allow_subgroup)

    invalid_reason = standardized_ig_symbol.get_reason_why_invalid(enforce_functional)

    if invalid_reason is None:
        return standardized_ig_symbol.compile(precision)

    if log_failures:
        _utils.warn_failure(
            reason_for_failure=invalid_reason,
            original_input=symbol,
            attempted_fix=standardized_ig_symbol.compile("allele"),
            species=species,
            logger=logger,
        )

    if on_fail == "reject":
        return None

    return symbol


def standardise(*args, **kwargs) -> Optional[str]:
    """
    Alias for :py:func:`tidytcells.ig.standardize`.
    """
    return standardize(*args, **kwargs)
