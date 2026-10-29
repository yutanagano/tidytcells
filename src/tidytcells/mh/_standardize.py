import logging
from tidytcells._config import config
from tidytcells import _utils
from tidytcells._utils import Parameter
from tidytcells._standardized_gene_symbol import (
    ImgtHlaSymbolStandardizer,
    MroHlaSymbolStandardizer,
    MroH2MhHaplotypeSymbolStandardizer,
    MroRT1MhHaplotypeSymbolStandardizer,
    MroBoLaSymbolStandardizer,
    MroMamuSymbolStandardizer,
    MroPatrSymbolStandardizer,
    MroSlaSymbolStandardizer
)
from typing import Dict, Optional, Type

from tidytcells.result._mh_gene import MhGene

logger = logging.getLogger(__name__)


SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS: Dict[str, Dict[str, Type]] = {
    "homosapiens": {"MRO": MroHlaSymbolStandardizer, "IPD": ImgtHlaSymbolStandardizer},
    "musmusculus": {"MRO": MroH2MhHaplotypeSymbolStandardizer},
    "macacamulatta": {"MRO": MroMamuSymbolStandardizer},
    "susscrofa": {"MRO": MroSlaSymbolStandardizer},
    "bostaurus": {"MRO": MroBoLaSymbolStandardizer},
    "pantroglodytes": {"MRO": MroPatrSymbolStandardizer},
    "rattusnorvegicus": {"MRO": MroRT1MhHaplotypeSymbolStandardizer},

}


def standardize(
    symbol: Optional[str] = None,
    species: Optional[str] = None,
    log_failures: Optional[bool] = None,
    gene: Optional[str] = None,
    suppress_warnings: Optional[bool] = None,
) -> MhGene:
    """
    Attempt to standardize an MH gene / allele symbol to be compliant with the names reported in the MHC Restriction Ontology (MRO).
    For most species, MRO follows IPD-MHC allele nomenclature (e.g. ``HLA-A*02:01``, ``Mamu-A1*001:01``).
    For *Mus musculus* and *Rattus norvegicus*, MRO names molecules by locus and haplotype (e.g. ``H2-Kb``, ``RT1-Aa``).

    MRO contains only functional MHC proteins. For ``"homosapiens"``, non-functional pseudogenes are also included when switching to the "IPD" database:
    tidytcells.set_mh_database("IPD")


    .. note::
        This function will only verify the validity of an MH gene/allele up to the level of the protein.
        Any further precise allele designations will not be verified, apart from the requirement that the format (colon-separated numbers) look valid.



    :param symbol:
        Potentially non-standardized MH gene / allele symbol.
    :type symbol:
        str
    :param species:
        Can be specified to standardize to a TR symbol that is known to be valid for that species (see above for supported species).
        Supported species are: ``"homosapiens"``, ``"musmusculus"``, ``"macacamulatta"``, ``"susscrofa"``, ``"bostaurus"``, ``"pantroglodytes"``, ``"rattusnorvegicus"`` or `"any"``
        Defaults to ``"homosapiens"``.
    :type species:
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


        Known alternative MRO names are mapped to the current MRO name, and mutant molecules are supported.

        >>> tt.mh.standardize("HLA-Cw*0301").symbol
        'HLA-C*03:04'
        >>> tt.mh.standardize("HLA-A0201 K66A mutant").symbol
        'HLA-A*02:01 K66A mutant'
        >>> tt.mh.standardize("HLA-A0201 K66A mutant").mutation
        'K66A'


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

            IF the configured database is not available for the species:
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
        tested_species = [species for species in SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS.keys()
                          if config.mhc_database in SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS[species].keys()]
        error = f'Failed with any of the following species: {', '.join(tested_species)}.'
        if config.mhc_database == "IPD":
            error += " To test more species, please switch to MRO using: tidytcells.set_mh_database('MRO')"
        best_attempt_result = MhGene(symbol, error)

        for species in tested_species:
            mh_standardizer = SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS[species][config.mhc_database](symbol)

            if mh_standardizer.result.is_standardized:
                return mh_standardizer.result

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

    if config.mhc_database not in SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS[species]:
        if log_failures:
            _utils.warn_unsupported_mhc_database(species, SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS, logger)
        return MhGene(symbol, f'Unsupported database: "{config.mhc_database}" for species "{species}"')

    standardizer_cls = SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS[species][config.mhc_database]
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