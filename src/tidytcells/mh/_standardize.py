import logging
from tidytcells import _utils
from tidytcells._utils import Parameter
from tidytcells._standardized_gene_symbol import (
    HlaSymbolStandardizer,
    MusMusculusMhSymbolStandardizer,
)
from typing import Dict, Optional, Type, Union

from tidytcells._utils.result import MhGene

logger = logging.getLogger(__name__)


SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS: Dict[str, Type[Union[HlaSymbolStandardizer, MusMusculusMhSymbolStandardizer]]] = {
    "homosapiens": HlaSymbolStandardizer,
    "musmusculus": MusMusculusMhSymbolStandardizer,
}


def standardize(
    symbol: Optional[str] = None,
    species: Optional[str] = None,
    database: Optional[str] = None,
    log_failures: Optional[bool] = None,
    gene: Optional[str] = None,
    suppress_warnings: Optional[bool] = None,
) -> MhGene:
    """
    Attempt to standardize an MH gene / allele symbol to be IMGT-compliant. # todo update docs once MRO mapping available

    .. topic:: Supported species

        - ``"homosapiens"``
        - ``"musmusculus"``

    .. note::
        This function will only verify the validity of an MH gene/allele up to the level of the protein.
        Any further precise allele designations will not be verified, apart from the requirement that the format (colon-separated numbers) look valid.
        The reasons for this is firstly because new alleles at that level are added to the IMGT list quite often and so accurate verification is difficult,
        secondly because people rarely need verification to such a precise level, and finally because such verification costs more computational effort with diminishing returns.

    :param symbol:
        Potentially non-standardized MH gene / allele symbol.
    :type symbol:
        str
    :param species:
        Can be specified to standardise to a TR symbol that is known to be valid for that species (see above for supported species).
        If set to ``"any"``, then first attempts standardisation for *Homo sapiens*, then *Mus musculus*.
        Defaults to ``"homosapiens"``.

        .. note::
            From version 3, the default behaviour will change to ``"any"``.

    :type species:
        str
    :param database:
        Which gene database to use. Defaults to ``"MRO"``, alternatively, ``"IMGT"`` can be selected.
        Note that IMGT uses a non-standard representation of mouse MH genes, and using MRO is therefore recommended.
        See also: https://github.com/IEDB/MRO
    :type database:
        str
    :param log_failures:
        Report standardisation failures through logging (at level ``WARNING``).
        Defaults to ``True``.
    :type log_failures:
        bool
    :param gene:
        Alias for `symbol`.
    :type gene:
        str
    :param suppress_warnings:
        Disable warnings that are usually logged when standardisation fails.
        Deprecated in favour of `log_failures`.
    :type suppress_warnings:
        bool

    :return:
        This method will return a MhGeneResult object with the following attributes:
            - success (bool): True if the standardiation was successful, False otherwise.
            - failed (bool): the inverse of success.
            - allele (str): the standardized symbol at the allele level, if standardisation to this level was successful, otherwise None.
            - protein (str): the standardized symbol at the protein level, if standardisation to this level was successful, otherwise None. This is only available for human data.
            - gene (str): the standardized symbol at the gene level, if standardisation was successful, otherwise None.
            - highest_precision (str): the most precise version of the standardized symbol (allele > gene) if standardisation was successful, otherwise None.
            - error (str): the error message, only if standardisation failed, otherwise None.
            - attempted_fix (str): the best attempt at fixing the input symbol, only of standardisation failed, otherwise None.
            - original_input (str): the original input symbol.
    :rtype:
        MhGeneResult


    .. topic:: Example usage

        MH standardised results will be returned as a MhGeneResult (or HLAGeneResult for human genes).
        When standardisation is a success, attributes 'allele', 'protein' and 'gene' can be used to retrieve the corrected information.

        >>> result = tt.mh.standardize("HLA-DRB3*01:01:02:01")
        >>> result.is_standardized
        True
        >>> result.allele
        'HLA-DRB3*01:01:02:01'
        >>> result.protein
        'HLA-DRB3*01:01'
        >>> result.gene
        'HLA-DRB3'


        Attributes 'allele', 'protein' and 'gene' only return a result if the symbol could be standardised up to that level.
        Attribute 'highest_precision' is never None for a successful standardisation, and always returns the most
        detailed available result between 'allele', 'protein' and 'gene'.

        >>> tt.mh.standardize("HLA-DRB3*01:01:02:01").symbol
        'HLA-DRB3*01:01:02:01'
        >>> tt.mh.standardize("HLA-DRB3").allele
        None
        >>> tt.mh.standardize("HLA-DRB3").gene
        'HLA-DRB3'
        >>> tt.mh.standardize("HLA-DRB3").symbol
        'HLA-DRB3'

        Non-standardised input strings will intelligently be corrected to IMGT-compliant symbols.

        >>> tt.mh.standardize("A1").allele
        'HLA-A*01'

        *Mus musculus* is a supported species. # todo update example with MRO

        >>> tt.mh.standardize("CRW2", species="musmusculus").gene
        'MH1-M5'

    .. topic:: Decision Logic

        #todo: update decision logic once mouse is updated with MRO mapping

        To provide an easy way to gauge the scope and limitations of standardisation, below is a simplified overview of the decision logic employed when attempting to standardize an MH symbol.
        For more detail, please refer to the `source code <https://github.com/yutanagano/tidytcells>`_.

        .. code-block:: none

            IF the specified species is not supported for standardisation:
                RETURN original symbol without modification

            ELSE:
                // attempt standardisation
                {
                    IF symbol is already in IMGT-compliant form:
                        set standardisation status as successful
                        skip rest of standardisation

                    IF symbol is a known deprecated symbol:
                        overwrite symbol with current IMGT-compliant symbol
                        set standardisation status as successful
                        skip rest of standardisation

                    // the rest is only applicable when species is set to homo sapiens
                    add "HLA-" to the beginning of the symbol if necessary                  //e.g. A -> HLA-A
                    replace "Cw" with "C"                                                   //e.g. HLA-Cw -> HLA-C
                    add back forgotten asterisks if necessary                               //e.g. HLA-A01 -> HLA-A*01
                    add back forgotten colons if necessary                                  //e.g. HLA-A*0101 -> HLA-A*01:01
                    If symbol is now in IMGT-compliant form:
                        set standardisation status as successful
                        skip rest of standardisation

                    try adding or subtracting leading zeros from allele designation numbers //e.g. HLA-A*001 -> HLA-A*01
                    If symbol is now in IMGT-compliant form:
                        set standardisation status as successful
                        skip rest of standardisation

                    set standardisation status as failed
                }

                RETURN MhGeneResult

    """
    symbol = (
        Parameter(symbol, "symbol")
        .resolve_with_alias(gene, "gene")
        .throw_error_if_not_of_type(str)
        .value
    )
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
    database = (
        Parameter(database, "species")
        .set_default("MRO")
        .throw_error_if_not_one_of("MRO", "IMGT")
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

    if species == "any":
        best_attempt_result = MhGene(symbol, f'Failed with any species')

        for (
            species,
            standardizer_cls,
        ) in SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS.items():
            mh_standardizer = standardizer_cls(symbol)

            if mh_standardizer.result.is_standardized:
                return mh_standardizer.result

            if species == "homosapiens":
                best_attempt_result = mh_standardizer.result

        if log_failures:
            _utils.warn_result_failure(
                result=best_attempt_result,
                logger=logger,
            )

        return best_attempt_result

    if species not in SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS:
        if log_failures:
            _utils.warn_unsupported_species(species, "MH", logger)
        return MhGene(symbol, f'Unsupported species: {species}')

    standardizer_cls = SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS[species]
    mh_standardizer = standardizer_cls(symbol)

    if (not mh_standardizer.result.is_standardized) and log_failures:
        _utils.warn_result_failure(
            result=mh_standardizer.result,
            logger=logger,
        )

    return mh_standardizer.result



def standardise(*args, **kwargs) -> Optional[str]:
    """
    Alias for :py:func:`tidytcells.mh.standardize`.
    """
    return standardize(*args, **kwargs)
