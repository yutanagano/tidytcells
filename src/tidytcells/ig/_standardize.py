import logging
from tidytcells import _utils
from tidytcells.result._receptor_gene import ReceptorGene
from tidytcells._utils import Parameter
from tidytcells._standardized_gene_symbol import (
    HomoSapiensIgSymbolStandardizer, MusMusculusIgSymbolStandardizer, ReceptorGeneSymbolStandardizer,
)
from typing import Dict, Optional, Type

logger = logging.getLogger(__name__)


SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS: Dict[str, Type[ReceptorGeneSymbolStandardizer]] = {
    "homosapiens": HomoSapiensIgSymbolStandardizer,
    "musmusculus": MusMusculusIgSymbolStandardizer,
}


def standardize(
    symbol: Optional[str] = None,
    species: Optional[str] = None,
    enforce_functional: Optional[bool] = None,
    allow_subgroup: Optional[bool] = None,
    log_failures: Optional[bool] = None,
    gene: Optional[str] = None,
    suppress_warnings: Optional[bool] = None,
) -> ReceptorGene:
    """
    Attempt to standardize a IG gene / allele symbol to be IMGT-compliant.

    .. topic:: Supported species

        - ``"homosapiens"``
        - ``"musmusculus"``

    :param symbol:
        Potentially non-standardized IG gene / allele symbol.
    :type symbol:
        str
    :param species:
        Can be specified to standardize to an IG symbol that is known to be valid for that species (see above for supported species).
        If set to ``"any"``, then first attempts standardization for *Homo sapiens*, then *Mus musculus*.
        Defaults to ``"homosapiens"``.

        .. note::
            From version 3, the default behaviour will change to ``"any"``.
    :type species:
        str
    :param enforce_functional:
        If ``True``, disallows IG genes / alleles that are recognised by IMGT but are marked as non-functional (ORF or pseudogene).
        Defaults to ``False``.
    :type enforce_functional:
        bool
    :param allow_subgroup:
        If ``True``, allows valid subgroups (as well as more specific gene/allele symbos) to pass standardization.
        If ``False``, the supplied symbol must point to at least a specific gene.
        Defaults to ``False``.
    :type allow_subgroup:
        bool
    :param log_failures:
        Report standardization failures through logging (at level ``WARNING``).
        Defaults to ``True``.
    :type log_failures:
        bool
    :param gene:
        Alias for the parameter `symbol`.
    :type gene:
        str
    :param suppress_warnings:
        Disable warnings that are usually logged when standardization fails.
        Deprecated in favour of `log_failures`.
    :type suppress_warnings:
        bool

    :return:
        A standardized receptor gene wrapped in a :py:class:`~tidytcells.result.ReceptorGene` object.
        For details on how to use this output, please refer to the class documentation.
    :rtype:
        `~tidytcells.result.ReceptorGene`

    .. topic:: Example usage

        IG standardized results will be returned as a :py:class:`~tidytcells.result.ReceptorGene`.
        When standardization is a success, attributes 'allele', 'gene' and 'subgroup' can be used to retrieve the corrected information.

        >>> result = tt.ig.standardize("IGHV1-2*01")
        >>> result.is_standardized
        True
        >>> result.allele
        'IGHV1-2*01'
        >>> result.gene
        'IGHV1-2'
        >>> result.subgroup
        'IGHV1'

        Attributes 'allele', 'gene' and 'subgroup' only return a result if the symbol could be standardized up to that level.
        Attribute 'symbol' is never None for a successful standardization, and always returns the most
        detailed available result between 'allele', 'gene' and 'subgroup'.

        >>> tt.ig.standardize("IGHV1-2*01").symbol
        'IGHV1-2*01'
        >>> tt.ig.standardize("IGHV1-2").allele
        None
        >>> tt.ig.standardize("IGHV1-2").gene
        'IGHV1-12'
        >>> tt.ig.standardize("IGHV1-2").symbol
        'IGHV1-12'

        Non-standardized input strings will intelligently be corrected to IMGT-compliant gene / allele symbols.

        >>> tt.ig.standardize("hj1").symbol
        'IGHJ1'

        The `enforce_functional` setting will cause non-functional genes or alleles to be rejected.
        For failed standardizations, the 'error' attribute explains why the standardization failed, and
        the 'attempted_fix' attribute contains the best attempted result found during standardization.

        >>> result = tt.ig.standardize("ighV1-12", enforce_functional=True)
        >>> result.is_standardized
        False
        >>> result.error
        'Gene has no functional alleles'
        >>> result.attempted_fix
        'IGHV1-12'

        Known synonyms are included in the standardization

        >>> tt.ig.standardize("A10").symbol
        'IGKV6D-21'

        *Mus musculus* is a supported species.

        >>> tt.ig.standardize("IGHV2-2", species="musmusculus").gene
        'IGHV2-2'

        Other available properties are 'original_input', 'species', 'receptor_type', 'locus' and 'gene_type'.

        >>> result = tt.ig.standardize("IGHV01-02")
        >>> result.symbol
        'IGHV1-2'
        >>> result.original_input
        'IGHV01-02'
        >>> result.species
        'homosapiens'
        >>> result.receptor_type
        'IG'
        >>> result.locus
        'IGH'
        >>> result.gene_type
        'V'

        Utility method 'get_all_alleles' can be used to retrieve all (functional) alleles for a given symbol.

        >>> result = tt.ig.standardize("IGHV1-3")
        >>> result.get_all_alleles()
        ['IGHV1-3*01', 'IGHV1-3*02', 'IGHV1-3*03', 'IGHV1-3*04', 'IGHV1-3*05']

        >>> result = tt.ig.standardize("IGHV1-67")
        >>> result.get_all_alleles(enforce_functional=True)
        []
        >>> result.get_all_alleles(enforce_functional=True)
        ['IGHV1-67*02', 'IGHV1-67*03', 'IGHV1-67*01']

        Utility method 'get_aa_sequences' can be used to retrieve known amino acid sequences per allele.
        Using sequence_type 'ALL' shows all available sequence data for each allele.

        >>> result = tt.ig.standardize("IGHV1-3*01")
        >>> result.get_aa_sequences(sequence_type="CDR1")
        {'IGHV1-3*01': 'GYTFTSYA'}
        >>> result.get_aa_sequences(sequence_type="CDR2")
        {'IGHV1-3*01': 'INAGNGNT'}
        >>> result = tt.ig.standardize("IGLJ3")
        >>> result.get_aa_sequences(sequence_type="ALL")
        {
        'IGLJ3*01': {
            'J-MOTIF': 'FGGG', 'J-REGION': 'VVFGGGTKLTVL', 'functionality': 'F'
            },
        'IGLJ3*02': {
            'J-MOTIF': 'FGGG', 'J-REGION': 'WVFGGGTKLTVL', 'functionality': 'F'
            }
        }


    .. topic:: Decision Logic

        To provide an easy way to gauge the scope and limitations of standardization, below is a simplified overview of the decision logic employed when attempting to standardize a TR symbol.
        For more detail, please refer to the `source code <https://github.com/yutanagano/tidytcells>`_.

        .. code-block:: none

            0. sanity-check input
            Skip standardization if invalid parameters are passed (invalid amino acids in sequence, invalid species, etc)

            1. attempt standardization
            IF symbol is already in IMGT-compliant form:
                set standardization status as successful, skip to step 2

            IF symbol is a known deprecated symbol:
                overwrite symbol with current IMGT-compliant symbol
                set standardization status as successful, skip to step 2

            replace "." with "-"                                        //e.g. IGHV1.2 -> IGHV1-2
            add back any missing backslashes                            //e.g. IGHV1OR15-1 -> IGHV1/OR15-1
            remove any unnecessary trailing zeros                       //e.g. IGHV1-02 -> IGHV1-2
            IF symbol is now in IMGT-compliant form:
                set standardization status as successful, skip to step 2

            add "IG" to the beginning of the symbol if necessary        //e.g. HV1-18 -> IGHV1-18
            IF symbol is now in IMGT-compliant form:
                set standardization status as successful, skip to step 2

            try removing "-1" from the end of the symbol                //e.g. IGHJ1-1 -> IGHJ1
            IF symbol is now in IMGT-compliant form:
                set standardization status as successful, skip to step 2

            set standardization status as failed

            2. finalisation
            IF standardization has not failed:
                consider standardization a success

            RETURN :py:class:`~tidytcells.result.ReceptorGene`
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
        best_attempt_result = ReceptorGene(symbol, f'Failed with any species')

        for (
            species,
            standardizer_cls,
        ) in SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS.items():
            ig_standardizer = standardizer_cls(symbol,
                                                        enforce_functional=enforce_functional,
                                                        allow_subgroup=allow_subgroup)

            if ig_standardizer.result.is_standardized:
                return ig_standardizer.result

            if species == "homosapiens":
                best_attempt_result = ig_standardizer.result

        if log_failures:
            _utils.warn_result_failure(
                result=best_attempt_result,
                logger=logger,
            )

        return best_attempt_result

    if species not in SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS:
        if log_failures:
            _utils.warn_unsupported_species(species, "IG", logger)
        return ReceptorGene(symbol, f'Unsupported species: {species}')

    standardizer_cls = SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS[species]
    ig_standardizer = standardizer_cls(symbol,
                                                enforce_functional=enforce_functional,
                                                allow_subgroup=allow_subgroup)

    if (not ig_standardizer.result.is_standardized) and log_failures:
        _utils.warn_result_failure(
            result=ig_standardizer.result,
            logger=logger,
        )

    return ig_standardizer.result



def standardise(*args, **kwargs) -> ReceptorGene:
    """
    Alias for :py:func:`tidytcells.ig.standardize`.

    :rtype:
        :py:class:`~tidytcells.result.ReceptorGene`
    """
    return standardize(*args, **kwargs)
