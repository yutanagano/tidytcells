import logging
from tidytcells import _utils
from tidytcells.result import ReceptorGene
from tidytcells._utils import Parameter
from tidytcells._standardized_gene_symbol import (
    ReceptorGeneSymbolStandardizer,
    HomoSapiensTrSymbolStandardizer,
    MusMusculusTrSymbolStandardizer,
)
from typing import Dict, Optional, Type

logger = logging.getLogger(__name__)


SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS: Dict[str, Type[ReceptorGeneSymbolStandardizer]] = {
    "homosapiens": HomoSapiensTrSymbolStandardizer,
    "musmusculus": MusMusculusTrSymbolStandardizer,
}


def standardize(
    symbol: Optional[str] = None,
    species: Optional[str] = None,
    enforce_functional: Optional[bool] = None,
    allow_subgroup: Optional[bool] = None,
    log_failures: Optional[str] = None,
    gene: Optional[str] = None,
    suppress_warnings: Optional[bool] = None,
) -> ReceptorGene:
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
        Can be specified to standardize to a TR symbol that is known to be valid for that species (see above for supported species).
        If set to ``"any"``, then first attempts standardization for *Homo sapiens*, then *Mus musculus*.
        Defaults to ``"homosapiens"``.

        .. note::
            From version 3, the default behaviour will change to ``"any"``.

    :type species:
        str
    :param enforce_functional:
        If ``True``, disallows TR genes / alleles that are recognised by IMGT but are marked as non-functional (ORF or pseudogene).
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
        :py:class:`~tidytcells.result.ReceptorGene`

    .. topic:: Example usage

        TR standardized results will be returned as a :py:class:`~tidytcells.result.ReceptorGene`.
        When standardization is a success, attributes 'allele', 'gene' and 'subgroup' can be used to retrieve the corrected information.

        >>> result = tt.tr.standardize("TRAV1-1*01")
        >>> result.is_standardized
        True
        >>> result.allele
        'TRAV1-1*01'
        >>> result.gene
        'TRAV1-1'
        >>> result.subgroup
        'TRAV1'

        Attributes 'allele', 'gene' and 'subgroup' only return a result if the symbol could be standardized up to that level.
        Attribute 'symbol' is never None for a successful standardization, and always returns the most
        detailed available result between 'allele', 'gene' and 'subgroup'.

        >>> tt.tr.standardize("TRAV1-1*01").symbol
        'TRAV1-1*01'
        >>> tt.tr.standardize("TRAV1-1").allele
        None
        >>> tt.tr.standardize("TRAV1-1").gene
        'TRAV1-1'
        >>> tt.tr.standardize("TRAV1-1").symbol
        'TRAV1-1'

        Non-standardized input strings will intelligently be corrected to IMGT-compliant gene / allele symbols.

        >>> tt.tr.standardize("aj1").gene
        'TRAJ1'

        The `enforce_functional` setting will cause non-functional genes or alleles to be rejected.
        For failed standardizations, the 'error' attribute explains why the standardization failed, and
        the 'attempted_fix' attribute contains the best attempted result found during standardization.

        >>> result = tt.tr.standardize("tcrBV1", enforce_functional=True)
        >>> result.is_standardized
        False
        >>> result.error
        'Gene has no functional alleles'
        >>> result.attempted_fix
        'TRBV1'

        Known synonyms are included in the standardization

        >>> tt.tr.standardize("V4P").symbol
        'TRGV11'

        *Mus musculus* is a supported species.

        >>> tt.tr.standardize("TRBV1", species="musmusculus").gene
        'TRBV1'

        Other available properties are 'original_input', 'species', 'receptor_type', 'locus' and 'gene_type'.

        >>> result = tt.tr.standardize("TRAV01-01")
        >>> result.symbol
        'TRAV1-1'
        >>> result.original_input
        'TRAV01-01'
        >>> result.species
        'homosapiens'
        >>> result.receptor_type
        'TR'
        >>> result.locus
        'TRA'
        >>> result.gene_type
        'V'

        Utility method 'get_all_alleles' can be used to retrieve all (functional) alleles for a given symbol.

        >>> result = tt.tr.standardize("TRAV1-1")
        >>> result.get_all_alleles()
        ['TRAV1-1*01', 'TRAV1-1*02']

        >>> result = tt.tr.standardize("TRAV15")
        >>> result.get_all_alleles(enforce_functional=True)
        []
        >>> result.get_all_alleles(enforce_functional=True)
        ['TRAV15*02', 'TRAV15*03', 'TRAV15*01']

        Utility method 'get_aa_sequences' can be used to retrieve known amino acid sequences per allele.
        Using sequence_type 'ALL' shows all available sequence data for each allele.

        >>> result = tt.tr.standardize("TRAV1-1*01")
        >>> result.get_aa_sequences(sequence_type="CDR1")
        {'TRAV1-1*01': 'TSGFYG'}
        >>> result.get_aa_sequences(sequence_type="CDR2")
        {'TRAV1-1*01': 'NALDGL'}
        >>> result = tt.tr.standardize("TRAJ15")
        >>> result.get_aa_sequences(sequence_type="ALL")
        {
        'TRAJ15*01': {
            'J-MOTIF': 'FGKG', 'J-REGION': 'NQAGTALIFGKGTTLSVSS', 'functionality': 'F'
            },
        'TRAJ15*02': {
            'J-MOTIF': 'FGKG', 'J-REGION': 'NQAGTALIFGKGTHLSVSS', 'functionality': 'F'
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
                set standardization status as successful, skip to step 2.

            replace "TCR" with "TR"                                     //e.g. TCRAV1-1 -> TRAV1-1
            replace "S" with "-"                                        //e.g. TRAV1S1 -> TRAV1-1
            replace "." with "-"                                        //e.g. TRAV1.1 -> TRAV1-1
            add back any missing backslashes                            //e.g. TRAV14DV4 -> TRAV14/DV4
            remove any unnecessary trailing zeros                       //e.g. TRAV1-01 -> TRAV1-1
            IF symbol is now in IMGT-compliant form:
                set standardization status as successful, skip to step 2


            add "TR" to the beginning of the symbol if necessary        //e.g. AV1-1 -> TRAV1-1
            IF symbol is now in IMGT-compliant form:
                set standardization status as successful, skip to step 2

            resolve compound TRAV/TRDV designation if necessary         //e.g. TRDV4 -> TRAV14/DV4 or TRAV14 -> TRAV14/DV4
            IF symbol is now in IMGT-compliant form:
                set standardization status as successful, skip to step 2

            try removing "-1" from the end of the symbol                //e.g. TRAV1-1 -> TRAV1
            IF symbol is now a valid IMGT-compliant *gene* (do not correct to subgroup):
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
            tr_standardizer = standardizer_cls(symbol,
                                                        enforce_functional=enforce_functional,
                                                        allow_subgroup=allow_subgroup)

            if tr_standardizer.result.is_standardized:
                return tr_standardizer.result

            if species == "homosapiens":
                best_attempt_result = tr_standardizer.result

        if log_failures:
            _utils.warn_result_failure(
                result=best_attempt_result,
                logger=logger,
            )

        return best_attempt_result

    if species not in SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS:
        if log_failures:
            _utils.warn_unsupported_species(species, "TR", logger)
        return ReceptorGene(symbol, f'Unsupported species: {species}')

    standardizer_cls = SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS[species]
    tr_standardizer = standardizer_cls(symbol,
                                                enforce_functional=enforce_functional,
                                                allow_subgroup=allow_subgroup)

    if (not tr_standardizer.result.is_standardized) and log_failures:
        _utils.warn_result_failure(
            result=tr_standardizer.result,
            logger=logger,
        )

    return tr_standardizer.result


def standardise(*args, **kwargs) -> ReceptorGene:
    """
    Alias for :py:func:`tidytcells.tr.standardize`.

    :rtype:
        :py:class:`~tidytcells.result.ReceptorGene`
    """
    return standardize(*args, **kwargs)
