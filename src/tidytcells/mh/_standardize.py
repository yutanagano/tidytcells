from typing import Dict, Optional, Type

from tidytcells import _utils
from tidytcells._utils import Parameter
from tidytcells._standardized_gene_symbol import (
    StandardizedGeneSymbol,
    StandardizedHlaSymbol,
    StandardizedMusMusculusMhSymbol,
)


SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS: Dict[str, Type[StandardizedGeneSymbol]] = {
    "homosapiens": StandardizedHlaSymbol,
    "musmusculus": StandardizedMusMusculusMhSymbol,
}


def standardize(
    gene: Optional[str] = None,
    species: str = "homosapiens",
    precision: str = "allele",
    on_fail: str = "reject",
    suppress_warnings: bool = False,
) -> tuple:
    """
    Attempt to standardize an MH gene name to be IMGT-compliant.

    .. topic:: Supported species

        - ``"homosapiens"``
        - ``"musmusculus"``

    .. note::
        This function will only verify the validity of an MH gene/allele up to the level of the protein.
        Any further precise allele designations will not be verified, apart from the requirement that the format (colon-separated numbers) look valid.
        The reasons for this is firstly because new alleles at that level are added to the IMGT list quite often and so accurate verification is difficult, secondly because people rarely need verification to such a precise level, and finally because such verification costs more computational effort with diminishing returns.

    :param gene:
        Potentially non-standardized MH gene name.
    :type gene:
        str
    :param species:
        Species to which the MH gene belongs (see above for supported species).
        Defaults to ``"homosapiens"``.
    :type species:
        str
    :param precision:
        The maximum level of precision to standardize to.
        ``"allele"`` standardizes to the maximum precision possible.
        ``"protein"`` keeps allele designators up to the level of the protein.
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
    :param suppress_warnings:
        Disable warnings that are usually emitted when standardisation fails.
        Defaults to ``False``.
    :type suppress_warnings:
        bool

    :return:
        If the specified ``species`` is supported, and ``gene`` could be standardized, then return the standardized gene name.
        If ``species`` is unsupported, then the function does not attempt to standardize, and returns the unaltered ``gene`` string.
        Else follows the behvaiour as set by ``on_fail``.
    :rtype:
        Union[str, None]

    .. topic:: Example usage

        Input strings will intelligently be corrected to IMGT-compliant gene symbols.

        >>> tt.mh.standardize("A1")
        'HLA-A*01'

        The ``precision`` setting can truncate unnecessary information.

        >>> tt.mh.standardize("HLA-A*01", precision="gene")
        'HLA-A'

        *Mus musculus* is a supported species.

        >>> tt.mh.standardize("CRW2", species="musmusculus")
        'MH1-M5'

    .. topic:: Decision Logic

        To provide an easy way to gauge the scope and limitations of standardization, below is a simplified overview of the decision logic employed when attempting to standardize an MH symbol.
        For more detail, please refer to the `source code <https://github.com/yutanagano/tidytcells>`_.

        .. code-block:: none

            IF the specified species is not supported for standardization:
                RETURN original gene symbol without modification

            ELSE:
                // attempt standardization
                {
                    IF gene symbol is already in IMGT-compliant form:
                        set standardization status as successful
                        skip rest of standardization

                    IF gene symbol is a known deprecated symbol:
                        overwrite gene symbol with current IMGT-compliant symbol
                        set standardization status as successful
                        skip rest of standardization

                    // the rest is only applicable when species is set to homo sapiens
                    add "HLA-" to the beginning of the gene symbol if necessary             //e.g. A -> HLA-A
                    replace "Cw" with "C"                                                   //e.g. HLA-Cw -> HLA-C
                    add back forgotten asterisks if necessary                               //e.g. HLA-A01 -> HLA-A*01
                    add back forgotten colons if necessary                                  //e.g. HLA-A*0101 -> HLA-A*01:01
                    If gene symbol is now in IMGT-compliant form:
                        set standardization status as successful
                        skip rest of standardization

                    try adding or subtracting leading zeros from allele designation numbers //e.g. HLA-A*001 -> HLA-A*01
                    If gene symbol is now in IMGT-compliant form:
                        set standardization status as successful
                        skip rest of standardization

                    set standardization status as failed
                }

                IF standardization status is set to successful:
                    RETURN standardized gene symbol

                ELSE:
                    IF on_fail is set to "reject":
                        RETURN None
                    IF on_fail is set to "keep":
                        RETURN original gene symbol without modification
    """
    Parameter(gene, "gene").throw_error_if_not_of_type(str)
    Parameter(species, "species").throw_error_if_not_of_type(str)
    Parameter(precision, "precision").throw_error_if_not_one_of(
        "allele", "protein", "gene"
    )
    Parameter(on_fail, "on_fail").throw_error_if_not_one_of("reject", "keep")
    Parameter(suppress_warnings, "suppress_warnings").throw_error_if_not_of_type(bool)

    species = _utils.clean_and_lowercase(species)

    species_is_supported = species in SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS
    if not species_is_supported:
        if not suppress_warnings:
            _utils.warn_unsupported_species(species, "MH")
        return gene

    StandardizedMhSymbolClass = SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS[species]
    standardized_mh_symbol = StandardizedMhSymbolClass(gene)

    invalid_reason = standardized_mh_symbol.get_reason_why_invalid()
    if invalid_reason is not None:
        if not suppress_warnings:
            _utils.warn_failure(
                reason_for_failure=invalid_reason,
                original_input=gene,
                attempted_fix=standardized_mh_symbol.compile("allele"),
                species=species,
            )
        if on_fail == "reject":
            return None
        return gene

    return standardized_mh_symbol.compile(precision)


def standardise(*args, **kwargs):
    """
    Alias for :py:func:`tidytcells.mh.standardize`.
    """
    return standardize(*args, **kwargs)
