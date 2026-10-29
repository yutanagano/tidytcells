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
    Attempt to standardize an MH gene / allele symbol to its name in the MHC Restriction Ontology (MRO).

    For most species, MRO follows IPD-MHC allele nomenclature (e.g. ``HLA-A*02:01``, ``Mamu-A1*001:01``).
    For *Mus musculus* and *Rattus norvegicus*, MRO names molecules by locus and haplotype (e.g. ``H2-Kb``, ``RT1-Aa``).

    .. topic:: Supported species and databases

        - ``"homosapiens"``: ``"MRO"`` (default), ``"IPD"``
        - ``"musmusculus"``: ``"MRO"``
        - ``"macacamulatta"``: ``"MRO"``
        - ``"susscrofa"``: ``"MRO"``
        - ``"bostaurus"``: ``"MRO"``
        - ``"pantroglodytes"``: ``"MRO"``
        - ``"rattusnorvegicus"``: ``"MRO"``

        The database is set for all MH standardization with ``tidytcells.set_mh_database()``.
        MRO covers the MH molecules that present antigens to T cells.
        For *Homo sapiens*, the ``"IPD"`` database (IPD-IMGT/HLA) additionally covers other genes in the HLA region,
        such as pseudogenes, HLA-DM, HLA-DO and MICA, as well as G and P groups.

    .. note::
        For species with numbered alleles (e.g. ``HLA-A*02:01``), this function only verifies an allele up to the level of the protein (the first two fields).
        Further fields are not verified, apart from the requirement that they look valid (numerical, at least two digits).
        For *Mus musculus* and *Rattus norvegicus*, symbols are verified as a whole (e.g. ``H2-Kb``).

    :param symbol:
        Potentially non-standardized MH gene / allele symbol.
    :type symbol:
        str
    :param species:
        Species to standardize the symbol for (see above for supported species).
        If set to ``"any"``, each species available for the configured database is tried in the order listed above,
        and the first successful result is returned.
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
        A standardized MH gene wrapped in a :py:class:`~tidytcells.result.MhGene` object
        For details on how to use this output, please refer to the class documentation.
    :rtype:
        :py:class:`~tidytcells.result.MhGene`

    .. topic:: Example usage

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
        >>> tt.mh.standardize("HLA-DRB3").symbol
        'HLA-DRB3'
        >>> tt.mh.standardize("HLA-DRB3").allele is None
        True

        Non-standardized input strings will intelligently be corrected to valid symbols.

        >>> tt.mh.standardize("A1").symbol
        'HLA-A*01'
        >>> tt.mh.standardize("hla-a*02:01").symbol
        'HLA-A*02:01'
        >>> tt.mh.standardize("D6S204").symbol
        'HLA-C'

        Known alternative names in MRO are mapped to the current name, and mutant molecules are supported.

        >>> tt.mh.standardize("HLA-Cw*0301").symbol
        'HLA-C*03:04'
        >>> result = tt.mh.standardize("HLA-A0201 K66A mutant")
        >>> result.symbol
        'HLA-A*02:01 K66A mutant'
        >>> result.allele
        'HLA-A*02:01'
        >>> result.mutation
        'K66A'

        Other species are supported, using their own nomenclature.

        >>> tt.mh.standardize("mamu-a1*00101", species="macacamulatta").symbol
        'Mamu-A1*001:01'
        >>> tt.mh.standardize("H-2Kb", species="musmusculus").symbol
        'H2-Kb'
        >>> tt.mh.standardize("I-Ab", species="musmusculus").symbol
        'H2-IAb'
        >>> tt.mh.standardize("rt1aa", species="rattusnorvegicus").symbol
        'RT1-Aa'

        With ``species="any"``, the species is detected from the symbol.

        >>> result = tt.mh.standardize("SLA-1*04:01", species="any")
        >>> result.symbol
        'SLA-1*04:01'
        >>> result.species
        'susscrofa'

        For failed standardizations, the 'error' attribute explains why the standardization failed, and
        the 'attempted_fix' attribute contains the best attempted result found during standardization.

        >>> result = tt.mh.standardize("HLA-DRB365")
        >>> result.is_standardized
        False
        >>> result.error
        'Nonexistent allele for recognized gene'
        >>> result.attempted_fix
        'HLA-DRB3*65'

        Other available properties are 'original_input' and 'species'.

        >>> result.original_input
        'HLA-DRB365'
        >>> result.species
        'homosapiens'

    .. topic:: Decision Logic

        To provide an easy way to gauge the scope and limitations of standardization, below is a simplified overview of the decision logic employed when attempting to standardize an MH symbol.
        For more detail, please refer to the `source code <https://github.com/yutanagano/tidytcells>`_.

        .. code-block:: none

            IF species is "any":
                try each species available for the configured database, in order
                RETURN the first successful result, or a failed result listing the species tried

            IF the specified species is not supported for standardization:
                RETURN original symbol and error message

            IF the configured database is not available for the species:
                RETURN original symbol and error message

            split off a trailing "<mutation> mutant" description if present         //e.g. HLA-A*02:01 K66A mutant

            // species with numbered alleles (e.g. *Homo sapiens*, *Macaca mulatta*)
            {
                IF symbol is already valid:
                    RETURN standardized symbol

                IF symbol is a known alternative name of an allele:
                    overwrite symbol with the current allele name                   //e.g. HLA-Cw*0301 -> HLA-C*03:04

                IF gene is a known deprecated gene name:
                    overwrite gene with the current gene name                       //e.g. D6S204 -> HLA-C

                add the species prefix if necessary                                 //e.g. A -> HLA-A
                replace "Cw" with "C" (*Homo sapiens* only)                         //e.g. HLA-Cw -> HLA-C
                add back forgotten asterisks if necessary                           //e.g. HLA-A01 -> HLA-A*01
                add back forgotten colons if necessary                              //e.g. HLA-A*0101 -> HLA-A*01:01
                try adding or removing leading zeros in the allele designation      //e.g. HLA-A*001 -> HLA-A*01
                try adding or removing a trailing "1" in the gene name              //e.g. HLA-DQA -> HLA-DQA1

                after each step, RETURN standardized symbol if it is now valid
            }

            // species named by locus and haplotype (mus musculus, rattus norvegicus)
            {
                IF symbol (or its known alternative name) is valid, ignoring case:
                    RETURN standardized symbol

                try matching while ignoring the species prefix and hyphens          //e.g. H-2Kb, I-Ab -> H2-Kb, H2-IAb
                    (case-sensitive if the symbol has no prefix)                    //e.g. Kb -> H2-Kb
            }

            IF no valid symbol was found:
                RETURN failed result with error message and attempted fix

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