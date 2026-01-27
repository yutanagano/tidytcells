import logging
from tidytcells import _utils
from tidytcells._resources import AMINO_ACIDS
from tidytcells._utils.result import Junction
from tidytcells._utils.alignment import get_is_valid_locus_gene_fn
from tidytcells._utils.parameter import Parameter
from tidytcells._standardized_junction import (
    JunctionStandardizer,
    HomoSapiensTrJunctionStandardizer,
    HomoSapiensIgJunctionStandardizer,
    MusMusculusTrJunctionStandardizer,
    MusMusculusIgJunctionStandardizer,
)

from typing import Dict, Optional, Type

logger = logging.getLogger(__name__)


SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS: Dict[str, Dict[str, Type[JunctionStandardizer]]] = {
    "homosapiens": {"TR": HomoSapiensTrJunctionStandardizer,
                    "IG": HomoSapiensIgJunctionStandardizer},
    "musmusculus": {"TR": MusMusculusTrJunctionStandardizer,
                    "IG": MusMusculusIgJunctionStandardizer},
}

def standardize(
    seq: str,
    locus: str,
    j_symbol: Optional[str] = None,
    v_symbol: Optional[str] = None,
    species: Optional[str] = None,
    allow_c_correction: Optional[bool] = None,
    allow_fw_correction: Optional[bool] = None,
    enforce_functional_v: Optional[bool] = None,
    enforce_functional_j: Optional[bool] = None,
    allow_v_reconstruction: Optional[bool] = None,
    allow_j_reconstruction: Optional[bool] = None,
    log_failures: Optional[bool] = None,
    suppress_warnings: Optional[bool] = None,
) -> Junction:
    """
    Corrects a given CDR3/Junction sequence into a valid
    and complete Junction sequence, based on alignment to V and J genes.
    This may include recovery of incorrectly trimmed amino acids,
    correction of sequencing errors in the conserved start and end positions,
    or trimming of unnecessary amino acids at the beginning and end of the sequence.

    :param seq:
        The junction sequence.
    :type seq:
        str
    :param locus:
        String value representing the locus (TRA, TRB, IGH, IGL, etc;
        TR or IG may be used if a more precise locus is unknown).
        This is used to select an applicable subset of V and J genes for junction correction.
    :type locus:
        str
    :param j_symbol:
        The TR/IG J symbol used to correct the end of the junction sequence.
        If a specific J allele is supplied (e.g., human TRAJ1*01),
        the junction is corrected according to the known sequence
        for this J allele. If a less precise gene or subgroup is given
        (e.g., human TRAJ23 which has multiple alleles),
        all associated allele sequences will be tested for the best alignment.
        If no J symbol is given, all J genes for the given species + locus will be tested.
    :type j_symbol:
        str
    :param v_symbol:
        The TR/IG V symbol used to correct the start of the junction sequence.
        If a specific V allele is supplied (e.g., human TRAV1-1*01),
        the junction is corrected according to the known sequence
        for this V allele. If a less precise gene or subgroup is given
        (e.g., human TRAV1-1 which has multiple alleles, or TRAV1 which has multiple genes),
        all associated allele sequences will be tested for the best alignment.
        If no V symbol is given, all V genes for the given species + locus will be tested.
    :type v_symbol:
        str
    :param species:
        The species that produced the underlying receptor. Defaults to ``homosapiens``.
    :type species:
        str
    :param allow_c_correction:
        Whether to allow the first amino acid in the input sequence to be corrected to 'C'
        if it is a potential sequencing error (only "W", "S", "R", "G", "Y", "F")
        and correction improves the V gene alignment. Defaults to ``False``.
    :type allow_c_correction:
        bool
    :param allow_fw_correction:
        Whether to allow the last amino acid in the input sequence to be corrected to 'F' or 'W'
        if it is a potential sequencing error (only "I", "L", "V", "Y", "S", "C", "G",  "R")
        and correction improves the J gene alignment. Defaults to ``False``.
    :type allow_fw_correction:
        bool
    :param enforce_functional_v:
        Only consider V genes which are annotated to be 'functional' (excluding ORFs and pseudogenes).
        Defaults to ``True``.
    :type enforce_functional_v:
        bool
    :param enforce_functional_j:
        Only consider V genes which are annotated to be 'functional' (excluding ORFs and pseudogenes).
        Due to its shorter size it can be difficult to determine whether a J gene is truly functional,
        and including ORFs/pseudogenes is therefore recommended for J genes but not for V genes.
        Defaults to ``False``.
    :type enforce_functional_j:
        bool
    :param allow_v_reconstruction:
        Whether to allow more than just the 1 conserved C amino acid to be re-constructed from the V gene.
        It is recommended to only set this value to True if V symbol information is supplied.
        Defaults to ``False``.
    :type allow_v_reconstruction:
        bool
    :param allow_j_reconstruction:
        Whether to allow more than just the 1 conserved F / W / C amino acid to be re-constructed from the J gene.
        It is recommended to only set this value to True if J symbol information is supplied.
        Defaults to ``False``.
    :type allow_j_reconstruction:
        bool
    :param log_failures:
        Report standardisation failures through logging (at level ``WARNING``).
        Defaults to ``True``.
    :type log_failures:
        bool
    :param suppress_warnings:
        Disable warnings that are usually logged when standardisation fails.
        Deprecated in favour of `log_failures`.
    :type suppress_warnings:
        bool

    :return:
        This method will return a JunctionResult object with the following attributes:
            - success (bool): True if the standardiation was successful, False otherwise.
            - failed (bool): the inverse of success.
            - junction (str): the IMGT-junction, including conserved leading C and trailing F / W / C if the standardiation was successful, otherwise None.
            - cdr3 (str): the IMGT-CDR3, excluding conserved leading C and trailing F / W / C if the standardiation was successful, otherwise None.
            - error (str): the error message, only if standardisation failed, otherwise None.
            - attempted_fix (str): the best attempt at fixing the input sequence, only of standardisation failed, otherwise None.
            - original_input (str): the original input sequence.
            - species (str): the species used for the gene lookup to validate the CDR3 junction.
    :rtype:
        JunctionResult

    .. topic:: Example usage

        Strings that look like junction sequences will be accepted and corrected to capitalised form.
        Use attributes 'junction' and 'cdr3' to retrieve the corrected sequences.

        >>> result = tt.junction.standardize("csadaf", locus="TR")
        >>> result.junction
        'CSADAF'
        >>> result.cdr3
        'SADA'
        >>> result.is_standardized
        True

        Strings that are valid amino acid sequences but do not start and end
        with the appropriate residues can be corrected based on V/J gene or locus information.
        When no gene information is provided, all possible genes for a given locus are tried.

        >>> tt.junction.standardize("sada", locus="TR").junction
        'CSADAF'

        When the provided sequence is too long, it will be trimmed

        >>> tt.junction.standardize("yicsadafg", locus="TR").junction
        'CSADAF'

        Sequences which cannot be standardized to junctions (no matches with V/J genes) will be rejected.
        The result object will provide an error message and attempted partial fix.

        >>> result = tt.junction.standardize("ASWEHGH", locus="TR")
        >>> print(result.junction)
        None
        >>> result.is_standardized
        False
        >>> result.error
        'J alignment unsuccessful; J side reconstruction unsuccessful.'
        >>> result.attempted_fix
        'CASWEHGH'

        The conserved trailing residue can be intelligently inferred if
        `j_symbol` is supplied.

        >>> tt.junction.standardize("CSADKLI", locus="TR" j_symbol="TRAJ38").junction
        'CSADKLIW'
        >>> tt.junction.standardize("CSADKLI", locus="TR" j_symbol="TRAJ37").junction
        'CSADKLIF'

        Missing conserved leading C residues can be inferred from the `v_symbol`, even
        when the provided sequence already starts with a C.

        >>> tt.junction.standardize("CSYAYVF", locus="IG" v_symbol="IGLV2-11").junction
        'CCSYAYVF'
        >>> tt.junction.standardize("CSYAYVF", locus="IG" v_symbol="IGLV2-11").cdr3
        'CSYAYV'
        >>> tt.junction.standardize("CSYAYVF", locus="IG" v_symbol="IGLV2-14").junction
        'CSYAYVF'

        Auotmatic correction of sequencing errors in the leading/trailing amino acid may be
        enabled with `allow_c_correction` and `allow_fw_correction`. This will only correct
        amino acids which could have occurred with 1 nucleotide difference in the codon.

        >>> tt.junction.standardize("WASSPGVFGANVLTF", locus="TR", allow_c_correction=True).junction
        'CASSPGVFGANVLTF'

        Reconstruction of more than 1 amino acid can be enabled with `allow_v_reconstruction`
        and `allow_j_reconstruction`. Since this option potentially adds many amino acids to the
        provided sequence, it is recommended to set `v_symbol` and `j_symbol` to as detailed
        information as possible to ensure correct results.

        >>> tt.junction.standardize("MRESENMD", locus="TR", \
                                    allow_v_reconstruction=True, allow_j_reconstruction=True, \
                                    v_symbol="TRAV14/DV4", j_symbol="TRAJ12").junction
        'CAMRESENMDSSYKLIF'


    .. topic:: Decision Logic

        To provide an easy way to gauge the scope and limitations of
        standardisation, below is a simplified overview of the decision logic
        employed when attempting to standardize a junction sequence. For more
        detail, please refer to the
        `source code <https://github.com/yutanagano/tidytcells>`_.

        .. code-block:: none


            0. sanity-check input
            Skip standardisation if invalid parameters are passed (invalid amino acids in sequence, invalid species, etc)

            1. select candidate V/J genes
            Retrieve all possible V and J genes based on the provided `locus`, `j_symbol`, `v_symbol` and `enforce_functional_v`, `enforce_functional_j`

            IF an allele-level V/J symbol is provided for which no sequence information is known:
                set standardisation status as failed (no sequence information known for <symbol>)

            2. align sequence to V/J genes
            Attempt to align to each of the retrieved V and J gene sequences through a sliding window approach:
                - Allow any number of amino acids to be removed from V gene C-terminus or J gene N-terminus (simulating junction site deletion)
                - Keep only V genes with at least 1 amino acid overlap with the sequence and with C as conserved amino acid
                - Keep only J genes with at least 2 amino acids overlap and at most 1 mismatch within the alignment
                - Keep only V and J alignments with the highest overlap score (+1 score for match, -1.5 penalty for mismatch, 0 score for terminal overhangs)

            IF  (`allow_c_correction` is True AND sequence start is in ("W", "S", "R", "G", "Y", "F")) OR
                (`allow_fw_correction` is True AND sequence end is in ("I", "L", "V", "Y", "S", "C", "G",  "R")):

                Repeat alignments with 'sequencing error corrected' version of the sequence (start -> C; end -> F, then end -> W)
                IF a higher alignment score is found with the corrected sequence compared to original:
                    overwrite the original sequence with the corrected sequence
                    keep only the best alignments associated with the corrected sequence

            IF no alignments are found:
                set standardisation status as failed (<V/J> alignment unsuccessful)

            3. correct J side
            FOR each successful J gene alignment:
                compute the corrected sequence according to this alignment
                IF the alignment matches perfectly up to the conserved amino acid:
                    no J-side correction needed, skip to step 4.

            IF any corrected sequences add only 1 terminal amino acid:
                keep only these, drop all other corrected sequences
            ELSE:
                'trimming' corrections are kept or corrections adding >1 amino acid if `allow_j_reconstruction` is True.

            IF `j_symbol` is None:
                keep only corrections ending with canonical F / W conserved amino acids

            IF one possible J side correction remains:
                keep only this correction, skip to step 4

            IF no J side corrections remain:
                set standardisation status as failed (J side reconstruction unsuccessful)

            IF multiple contradicting J side corrections remain:
                set standardisation status as failed (J side reconstruction ambiguous)

            4. correct V side
            FOR each successful V gene alignment:
                compute the corrected sequence according to this alignment
                IF any alignment matches perfectly up to the conserved leading C:
                    no V-side correction needed, skip to step 4.

            IF any alignment adds only 1 leading C:
                keep only these, drop all other corrected sequences
            ELSE:
                'trimming' corrections are kept or corrections adding >1 amino acid if `allow_c_reconstruction` is True.

            IF one possible V side correction remains:
                keep only this correction, skip to step 5

            IF no V side corrections remain,:
                set standardisation status as failed (V side reconstruction unsuccessful)

            IF multiple contradicting V side corrections remain:
                set standardisation status as failed (V side reconstruction ambiguous)

            5. check length
            IF the length of the remaining sequence is < 6:
                set standardisation status as failed (junction too short)

            6. finalisation
            IF standardisation has not failed:
                consider standardisation a success
            RETURN JunctionResult

    """
    seq = Parameter(seq, "seq").throw_error_if_not_of_type(str).value
    locus = (
        Parameter(locus, "locus")
        .throw_error_if_not_one_of("TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL", "IG", "TR")
        .value
    )
    j_symbol = (
        Parameter(j_symbol, "j_symbol")
        .throw_error_if_not_of_type(str, optional=True)
        .throw_error_if_failed_test(test=get_is_valid_locus_gene_fn(locus, "J"),
                                    mssg=f"is not a valid J gene for \"locus\" {locus}",
                                    optional=True)
        .value
    )
    v_symbol = (
        Parameter(v_symbol, "v_symbol")
        .throw_error_if_not_of_type(str, optional=True)
        .throw_error_if_failed_test(test=get_is_valid_locus_gene_fn(locus, "V"),
                                    mssg=f"is not a valid V gene for \"locus\" {locus}",
                                    optional=True)
        .value
    )
    species = (
        Parameter(species, "species")
        .set_default("homosapiens")
        .throw_error_if_not_of_type(str)
        .value
    )
    allow_c_correction = (
        Parameter(allow_c_correction, "allow_c_correction")
        .set_default(False)
        .throw_error_if_not_of_type(bool)
        .value
    )
    allow_fw_correction = (
        Parameter(allow_fw_correction, "allow_fw_correction")
        .set_default(False)
        .throw_error_if_not_of_type(bool)
        .value
    )
    enforce_functional_v = (
        Parameter(enforce_functional_v, "enforce_functional_v")
        .set_default(True)
        .throw_error_if_not_of_type(bool)
        .value
    )
    enforce_functional_j = (
        Parameter(enforce_functional_j, "enforce_functional_j")
        .set_default(False)
        .throw_error_if_not_of_type(bool)
        .value
    )
    allow_v_reconstruction = (
        Parameter(allow_v_reconstruction, "allow_v_reconstruction")
        .set_default(False)
        .throw_error_if_not_of_type(bool)
        .value
    )
    allow_j_reconstruction = (
        Parameter(allow_j_reconstruction, "allow_j_reconstruction")
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
    original_input = seq

    seq = seq.upper().strip()

    for char in seq:
        if char not in AMINO_ACIDS:
            if log_failures:
                logger.warning(
                    f'Failed to standardize {original_input}. Not a valid amino acid sequence, found: {char}'
                )
            return Junction(original_input, f'Not a valid amino acid sequence, found: {char}')


    if species not in SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS:
        if log_failures:
            _utils.warn_unsupported_species(species, "junction", logger)

        return Junction(original_input, f'Unsupported species: {species}')


    if locus[0:2] not in SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS[species]:
        if log_failures:
            logger.warning(
                f'Unsupported locus: "{locus}" for species "{species}". ' f"Skipping {type} standardisation."
            )

        return Junction(original_input, f'Unsupported locus: "{locus}" for species "{species}"')


    standardizer_cls = SUPPORTED_SPECIES_AND_THEIR_STANDARDIZERS[species][locus[0:2]]
    result = standardizer_cls(seq=seq, locus=locus, j_symbol=j_symbol, v_symbol=v_symbol,
                                                      allow_c_correction=allow_c_correction,
                                                      allow_fw_correction=allow_fw_correction,
                                                      enforce_functional_v=enforce_functional_v,
                                                      enforce_functional_j=enforce_functional_j,
                                                      allow_v_reconstruction=allow_v_reconstruction,
                                                      allow_j_reconstruction=allow_j_reconstruction).result

    if (not result.is_standardized) and log_failures:
        _utils.warn_result_failure(
            result=result,
            logger=logger,
        )

    return result


def standardise(*args, **kwargs) -> Optional[str]:
    """
    Alias for :py:func:`tidytcells.junction.standardize`.
    """
    return standardize(*args, **kwargs)
