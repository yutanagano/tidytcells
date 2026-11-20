import logging
from tidytcells import config
from tidytcells import _utils
from tidytcells._standardized_gene_symbol import IgSymbolStandardizer
from tidytcells.result._receptor_gene import ReceptorGene
from tidytcells._utils import Parameter
from tidytcells._resources import IG_SYMBOLS_BY_SP_DB

logger = logging.getLogger(__name__)


def standardize(
    symbol: str | None = None,
    species: str | None = None,
    enforce_functional: bool | None = None,
    allow_subgroup: bool | None = None,
    log_failures: bool | None = None,
    gene: str | None = None,
    suppress_warnings: bool | None = None,
) -> ReceptorGene:
    """
    Attempt to standardize a IG gene / allele symbol to be IMGT-compliant.

    .. topic:: Supported species

        - ``"homosapiens"``
        - ``"musmusculus"``

    Parameters
    ----------
    symbol
        Potentially non-standardized IG gene / allele symbol.
    species
        Can be specified to standardize to an IG symbol that is known to be valid for that species (see above for supported species).
        Defaults to ``"homosapiens"``.
    enforce_functional
        If ``True``, disallows IG genes / alleles that are recognised by IMGT but are marked as non-functional (ORF or pseudogene).
        Defaults to ``False``.
    allow_subgroup
        If ``True``, allows valid subgroups (as well as more specific gene/allele symbos) to pass standardization.
        If ``False``, the supplied symbol must point to at least a specific gene.
        Defaults to ``False``.
    log_failures
        Report standardization failures through logging (at level ``WARNING``).
        Defaults to ``True``.
    gene
        Alias for the parameter `symbol`.
    suppress_warnings
        Disable warnings that are usually logged when standardization fails.
        Deprecated in favour of `log_failures`.

    Returns
    -------
    ~tidytcells.result.ReceptorGene
        A standardized receptor gene.
        For details on how to use this output, please refer to the class documentation.

    Examples
    --------

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
        >>> result.get_all_alleles(enforce_functional=False)
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
            'functionality': 'F', 'J-REGION': 'VVFGGGTKLTVL', 'J-MOTIF': 'FGGG'
            },
        'IGLJ3*02': {
            'functionality': 'F', 'J-REGION': 'WVFGGGTKLTVL', 'J-MOTIF': 'FGGG'
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

    if species not in IG_SYMBOLS_BY_SP_DB:
        msg = _utils.unsupported_species_message(species, IG_SYMBOLS_BY_SP_DB.keys())
        if log_failures:
            logger.warning(f"{msg} Skipping IG standardization.")
        return ReceptorGene(symbol, msg)

    receptor_databases = IG_SYMBOLS_BY_SP_DB[species]

    receptor_database = config.get_receptor_database()
    if receptor_database not in receptor_databases:
        msg = _utils.unsupported_receptor_database_message(
            receptor_database,
            species,
            receptor_databases.keys(),
        )
        if log_failures:
            logger.warning(f"{msg} Skipping IG standardization.")
        return ReceptorGene(symbol, msg)

    valid_gene_dictionary = receptor_databases[receptor_database]["symbols"]
    synonym_dictionary = receptor_databases[receptor_database].get("synonyms", {})
    germline_database_info = receptor_databases[receptor_database]["provenance"]

    ig_standardizer = IgSymbolStandardizer(
        symbol,
        enforce_functional=enforce_functional,
        allow_subgroup=allow_subgroup,
        species=species,
        germline_database_info=germline_database_info,
        valid_gene_dictionary=valid_gene_dictionary,
        synonym_dictionary=synonym_dictionary,
    )

    if (not ig_standardizer.result.is_standardized) and log_failures:
        logger.warning(_utils.standardization_failure_message(ig_standardizer.result))

    return ig_standardizer.result


def standardise(*args, **kwargs) -> ReceptorGene:
    """Alias for :py:func:`tidytcells.ig.standardize`."""
    return standardize(*args, **kwargs)
