import logging
from tidytcells import _utils
from tidytcells._utils import Parameter
from tidytcells._standardized_gene_symbol import (
    ImgtHlaSymbolStandardizer,
    MroHlaSymbolStandardizer,
    ImgtMusMusculusMhSymbolStandardizer,
    MroMusMusculusMhSymbolStandardizer,
)
from typing import Dict, Optional, Type

from tidytcells.result._mh_gene import MhGene

logger = logging.getLogger(__name__)


SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS: Dict[str, Dict[str, Type]] = {
    "homosapiens": {"IMGT": ImgtHlaSymbolStandardizer, "MRO": MroHlaSymbolStandardizer},
    "musmusculus": {"IMGT": ImgtMusMusculusMhSymbolStandardizer, "MRO": MroMusMusculusMhSymbolStandardizer},
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
    Attempt to standardize an MH gene / allele symbol to be MRO- or IMGT-compliant.

    .. topic:: Supported species and databases

        - ``"homosapiens"``: ``"MRO"``, ``"IMGT"``
        - ``"musmusculus"``: ``"MRO"``, ``"IMGT"``

    .. note::
        This function will only verify the validity of an MH gene/allele up to the level of the protein.
        Any further precise allele designations will not be verified, apart from the requirement that the format (colon-separated numbers) look valid.
        The reasons for this is firstly because new alleles at that level are added quite often and so accurate verification is difficult,
        secondly because people rarely need verification to such a precise level, and finally because such verification costs more computational effort with diminishing returns.

    :param symbol:
        Potentially non-standardized MH gene / allele symbol.
    :type symbol:
        str
    :param species:
        Can be specified to standardize to a TR symbol that is known to be valid for that species (see above for supported species).
        If set to ``"any"``, then first attempts standardization for *Homo sapiens*, then *Mus musculus*.
        Defaults to ``"homosapiens"``.

        .. note::
            From version 3, the default behaviour will change to ``"any"``.

    :type species:
        str
    :param database:
        Which database to standardize against: ``"MRO"`` (https://github.com/IEDB/MRO) or ``"IMGT"``.
        If species is ``"any"``, only species for which the database is available are attempted.
        Note that IMGT uses a non-standard representation of mouse MH genes (e.g. ``MH1-K1`` rather than ``H2-Kb``).
        Defaults to ``"MRO"``.
    :type database:
        str
    :param log_failures:
        Report standardization failures through logging (at level ``WARNING``).
        Defaults to ``True``.
    :type log_failures:
        bool
    :param gene:
        Alias for `symbol`.
    :type gene:
        str
    :param suppress_warnings:
        Disable warnings that are usually logged when standardization fails.
        Deprecated in favour of `log_failures`.
    :type suppress_warnings:
        bool

    :return:
        A standardized MHC gene wrapped in a :py:class:`~tidytcells.result.MhGene` object.
        For details on how to use this output, please refer to the class documentation.
    :rtype:
        :py:class:`~tidytcells.result.MhGene`


    .. topic:: Example usage

        MH standardized results will be returned as a :py:class:`~tidytcells.result.MhGene` (or :py:class:`~tidytcells.result.HLAGene` for human genes).
        When standardization is a success, attributes 'allele', 'protein' and 'gene' can be used to retrieve the corrected information.

        >>> result = tt.mh.standardize("HLA-DRB3*01:01:02:01")
        >>> result.is_standardized
        True
        >>> result.allele
        'HLA-DRB3*01:01:02:01'
        >>> result.protein
        'HLA-DRB3*01:01'
        >>> result.gene
        'HLA-DRB3'


        Attributes 'allele', 'protein' and 'gene' only return a result if the symbol could be standardized up to that level.
        Attribute 'symbol' is never None for a successful standardization, and always returns the most
        detailed available result between 'allele', 'protein' and 'gene'.

        >>> tt.mh.standardize("HLA-DRB3*01:01:02:01").symbol
        'HLA-DRB3*01:01:02:01'
        >>> tt.mh.standardize("HLA-DRB3").allele
        None
        >>> tt.mh.standardize("HLA-DRB3").gene
        'HLA-DRB3'
        >>> tt.mh.standardize("HLA-DRB3").symbol
        'HLA-DRB3'

        Non-standardized input strings will intelligently be corrected to compliant symbols.

        >>> tt.mh.standardize("A1").allele
        'HLA-A*01'

        *Mus musculus* is a supported species.

        >>> tt.mh.standardize("H-2Kb", species="musmusculus").symbol
        'H2-Kb'
        >>> tt.mh.standardize("I-Ab", species="musmusculus").symbol
        'H2-IAb'
        >>> tt.mh.standardize("CRW2", species="musmusculus", database="IMGT").gene
        'MH1-M5'

        Known alternative MRO names are mapped to the current MRO name, and mutant molecules are supported.

        >>> tt.mh.standardize("HLA-Cw*0301").symbol
        'HLA-C*03:04'
        >>> tt.mh.standardize("HLA-A0201 K66A mutant").symbol
        'HLA-A*02:01 K66A mutant'
        >>> tt.mh.standardize("HLA-A0201 K66A mutant").mutation
        'K66A'

        IMGT can be selected instead of MRO, which may give a different result for the same input.

        >>> tt.mh.standardize("I-Ab", species="musmusculus").symbol
        'H2-IAb'
        >>> tt.mh.standardize("I-Ab", species="musmusculus", database="IMGT").symbol
        'MH2-AB'


        For failed standardizations, the 'error' attribute explains why the standardization failed, and
        the 'attempted_fix' attribute contains the best attempted result found during standardization.

        >>> result = tt.mh.standardize(symbol="HLA-DRB365")
        >>> result.is_standardized
        False
        >>> result.error
        'Nonexistent allele for recognized gene'
        >>> result.attempted_fix
        'HLA-DRB3*65'

        Other available properties are 'original_input', 'species'.
        >>> result.original_input
        'HLA-DRB365'
        >>> result.original_input
        'homosapiens'

    .. topic:: Decision Logic

        To provide an easy way to gauge the scope and limitations of standardization, below is a simplified overview of the decision logic employed when attempting to standardize an MH symbol.
        For more detail, please refer to the `source code <https://github.com/yutanagano/tidytcells>`_.

        .. code-block:: none

            IF the specified species is not supported for standardization:
                RETURN original symbol and error message

            IF the specified database is not available for the species:
                RETURN original symbol and error message

            ELSE:
                // attempt standardization
                {
                    split off a trailing "<mutation> mutant" description if present       //e.g. HLA-A*02:01 K66A mutant

                IF symbol is a known alternative name of an allele:
                    overwrite symbol with the current allele name                     //e.g. HLA-Cw*0301 -> HLA-C*03:04

                // for mus musculus, where the database names haplotypes
                try matching while ignoring an H2 / H-2 prefix and hyphens              //e.g. H-2Kb, I-Ab -> H2-Kb, H2-IAb

                    IF symbol is already in compliant form:
                        set standardization status as successful
                        skip rest of standardization

                    IF database is IMGT and gene is a known deprecated symbol:
                        overwrite gene with current IMGT-compliant symbol
                        set standardization status as successful
                        skip rest of standardization

                    // the rest is only applicable when species is set to homo sapiens
                    add "HLA-" to the beginning of the symbol if necessary                  //e.g. A -> HLA-A
                    replace "Cw" with "C"                                                   //e.g. HLA-Cw -> HLA-C
                    add back forgotten asterisks if necessary                               //e.g. HLA-A01 -> HLA-A*01
                    add back forgotten colons if necessary                                  //e.g. HLA-A*0101 -> HLA-A*01:01
                    If symbol is now in compliant form:
                        set standardization status as successful
                        skip rest of standardization

                    try adding or subtracting leading zeros from allele designation numbers //e.g. HLA-A*001 -> HLA-A*01
                    If symbol is now in compliant form:
                        set standardization status as successful
                        skip rest of standardization

                    set standardization status as failed
                }

                RETURN :py:class:`~tidytcells.result.MhGene`

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
    database = (
        Parameter(database, "database")
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
            standardizers,
        ) in SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS.items():
            if database not in standardizers:
                continue

            mh_standardizer = standardizers[database](symbol)

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

    if database not in SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS[species]:
        if log_failures:
            logger.warning(
                f'Unsupported database: "{database}" for species "{species}" (valid options are: {list(SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS[species])}). Skipping MH standardization.'
            )
        return MhGene(symbol, f'Unsupported database: "{database}" for species "{species}"')

    standardizer_cls = SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS[species][database]
    mh_standardizer = standardizer_cls(symbol)

    if (not mh_standardizer.result.is_standardized) and log_failures:
        _utils.warn_result_failure(
            result=mh_standardizer.result,
            logger=logger,
        )

    return mh_standardizer.result



def standardise(*args, **kwargs) -> MhGene:
    """
    Alias for :py:func:`tidytcells.mh.standardize`.

    :rtype:
        :py:class:`~tidytcells.result.MhGene`
    """
    return standardize(*args, **kwargs)