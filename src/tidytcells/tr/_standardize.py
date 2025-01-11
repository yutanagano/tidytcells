import logging
from tidytcells import _utils
from tidytcells._utils import Parameter
from tidytcells._standardized_gene_symbol import (
    StandardizedSymbol,
    StandardizedHomoSapiensTrSymbol,
    StandardizedMusMusculusTrSymbol,
)
from typing import Dict, Optional, Type, Literal


logger = logging.getLogger(__name__)


SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS: Dict[str, Type[StandardizedSymbol]] = {
    "homosapiens": StandardizedHomoSapiensTrSymbol,
    "musmusculus": StandardizedMusMusculusTrSymbol,
}


def standardize(
    symbol: Optional[str] = None,
    species: Optional[str] = None,
    enforce_functional: Optional[bool] = None,
    precision: Optional[Literal["allele", "gene"]] = None,
    on_fail: Optional[Literal["reject", "keep"]] = None,
    log_failures: Optional[str] = None,
    gene: Optional[str] = None,
    suppress_warnings: Optional[bool] = None,
) -> str:
    """
    Attempt to standardize a TR gene / allele symbol to be IMGT-compliant.

    .. topic:: Supported species

        - ``"homosapiens"``
        - ``"musmusculus"``

    :param symbol:
        Potentially non-standardized TR gene / allele symbol.
    :type symbol:
        str
    :param species:
        Species to which the TR gene / allele belongs (see above for supported species).
        Defaults to ``"homosapiens"``.
    :type species:
        str
    :param enforce_functional:
        If ``True``, disallows TR genes / alleles that are recognised by IMGT but are marked as non-functional (ORF or pseudogene).
        Defaults to ``False``.
    :type enforce_functional:
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
        Union[str, None]

    .. topic:: Example usage

        Input strings will intelligently be corrected to IMGT-compliant gene / allele symbols.

        >>> tt.tr.standardize("aj1")
        'TRAJ1'

        The `precision` setting can truncate unnecessary information.

        >>> tt.tr.standardize("TRBV6-4*01", precision="gene")
        'TRBV6-4'

        The `enforce_functional` setting will cause non-functional genes or alleles to be rejected.

        >>> result = tt.tr.standardize("TRBV1", enforce_functional=True)
        Failed to standardize "TRBV1" for species homosapiens: gene has no functional alleles. Attempted fix "TRBV1".
        >>> print(result)
        None

        *Mus musculus* is a supported species.

        >>> tt.tr.standardize("TCRBV22S1A2N1T", species="musmusculus")
        'TRBV2'

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

                    replace "TCR" with "TR"                                     //e.g. TCRAV1-1 -> TRAV1-1
                    replace "S" with "-"                                        //e.g. TRAV1S1 -> TRAV1-1
                    replace "." with "-"                                        //e.g. TRAV1.1 -> TRAV1-1
                    add back any missing backslashes                            //e.g. TRAV14DV4 -> TRAV14/DV4
                    remove any unnecessary trailing zeros                       //e.g. TRAV1-01 -> TRAV1-1
                    IF symbol is now in IMGT-compliant form:
                        set standardization status as successful
                        skip rest of standardization


                    add "TR" to the beginning of the symbol if necessary   //e.g. AV1-1 -> TRAV1-1
                    IF symbol is now in IMGT-compliant form:
                        set standardization status as successful
                        skip rest of standardization

                    resolve compound TRAV/TRDV designation if necessary         //e.g. TRDV4 -> TRAV14/DV4 or TRAV14 -> TRAV14/DV4
                    IF symbol is now in IMGT-compliant form:
                        set standardization as successful
                        skip rest of standardization

                    try adding or removing "-1" from the end of the symbol //e.g. TRAV1 -> TRAV1-1
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
    precision = (
        Parameter(precision, "precision")
        .set_default("allele")
        .throw_error_if_not_one_of("allele", "gene")
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

    species = _utils.clean_and_lowercase(species)

    species_is_supported = species in SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS
    if not species_is_supported:
        if log_failures:
            _utils.warn_unsupported_species(species, "TR", logger)
        return symbol

    StandardizedTrSymbolClass = SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS[species]
    standardized_tr_symbol = StandardizedTrSymbolClass(symbol)

    invalid_reason = standardized_tr_symbol.get_reason_why_invalid(enforce_functional)
    if invalid_reason is not None:
        if log_failures:
            _utils.warn_failure(
                reason_for_failure=invalid_reason,
                original_input=symbol,
                attempted_fix=standardized_tr_symbol.compile("allele"),
                species=species,
                logger=logger,
            )
        if on_fail == "reject":
            return None
        return symbol

    return standardized_tr_symbol.compile(precision)


def standardise(*args, **kwargs):
    """
    Alias for :py:func:`tidytcells.tr.standardize`.
    """
    return standardize(*args, **kwargs)
