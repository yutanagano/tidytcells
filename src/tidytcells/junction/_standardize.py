import logging
import re
from tidytcells import aa, _utils
from typing import Literal, Optional
from tidytcells._utils.parameter import Parameter
from tidytcells._utils.conserved_aa_lookup import get_conserved_aa
from tidytcells._utils.trimming import process_junction

logger = logging.getLogger(__name__)



def standardize(
    seq: str,
    locus: Optional[str] = None,
    j_symbol: Optional[str] = None,
    species: Optional[str] = None,
    strict: Optional[bool] = None,
    j_strict: Optional[bool] = None,
    trimming: Optional[bool] = None,
    on_fail: Optional[Literal["reject", "keep"]] = None,
    log_failures: Optional[bool] = None,
    suppress_warnings: Optional[bool] = None,
) -> Optional[str]:
    """
    Ensures that a string value looks like a valid junction (CDR3) amino acid sequence.
    This function is a special variant of :py:func:`tidytcells.aa.standardize`.

    A valid junction sequence must:

    1. Be a valid amino acid sequence
    2. Begin with a cysteine (C)
    3. End with a phenylalanine (F) or a tryptophan (W) if no J symbol is supplied;
        or non-canonical cysteine (C) only if the supplied j_symbol is TRAJ35*01 and species is human

    :param seq:
        String value representing a junction sequence.
    :type seq:
        str
    :param locus:
        String value representing the locus (TRA, TRB, IGH, IGL, etc). If supplied, this is used to inform conserved trailing
        amino acid identification and trimming.
    :type locus:
        str
    :param j_symbol:
        String value representing the J symbol used to determine the correct conserved trailing amino acid sequence (F or W).
        This is determined based on the allele sequence. If less-precise gene or subgroup information is provided, an attempt will be made
        to extend the J symbol name to all known alleles.
        If not provided, the conserved trailing amino acid which may be added to the end of the sequence will default to 'F'.
    :type j_symbol:
        str
    :param species:
        String value representing the species of the J symbol. Defaults to ``homosapiens``.
    :type species:
        str
    :param strict:
        If ``True``, any string that does not look like a junction sequence is rejected.
        If ``False``, any inputs that are valid amino acid sequences but do not start with C and end with F/W are not rejected and instead are corrected by having a C appended to the beginning and an F appended at the end.
        Defaults to ``False``.
    :type strict:
        bool
    :param j_strict:
        If ``True``, the input will be rejected if a valid conserved trailing amino acid (F/W) cannot with certainty be determined from the J symbol.
        If ``False``, the default trailing amino acid F will be added to the end of the sequence if it cannot be determined from the J symbol.
        Defaults to ``False``.
    :type j_strict:
        bool
    :param trimming:
        If ``True``, the start and end of the supplied sequence may be trimmed until the conserved leading C and trailing F/W(/C).
        If ``False``, no trimming will be done.
        Defaults to ``False``.
    :type trimming:
        bool
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
    :param suppress_warnings:
        Disable warnings that are usually logged when standardisation fails.
        Deprecated in favour of `log_failures`.
    :type suppress_warnings:
        bool

    :return:
        If possible, a standardized version of the input string is returned.
        If the input string cannot be standardized, the function follows the behaviour as set by `on_fail`.
    :rtype:
        Optional[str]

    .. topic:: Example usage

        Strings that look like junction sequences will be accepted, and returned in capitalised form.

        >>> tt.junction.standardize("csadaf")
        'CSADAF'

        Strings that are valid amino acid sequences but do not stard and end with the appropriate residues will have a C and an F appended to its beginning and end as required.

        >>> tt.junction.standardize("sada")
        'CSADAF'

        However, setting `strict` to ``True`` will cause these cases to be rejected.

        >>> result = tt.junction.standardize("sada", strict=True)
        Input sadaf was rejected as it is not a valid junction sequence.
        >>> print(result)
        None

        By default, the conserved trailing amino acid is presumed to be 'F'.
        When `j_symbol` is provided, the correct trailing amino acid (F or W) is determined based on the symbol.

        >>> tt.junction.standardize("AELNAGNNRKLI")
        'CAELNAGNNRKLIF'

        >>> tt.junction.standardize("AELNAGNNRKLI", j_symbol="TRAJ38*01")
        'CAELNAGNNRKLIW'


    .. topic:: Decision Logic

        To provide an easy way to gauge the scope and limitations of standardization, below is a simplified overview of the decision logic employed when attempting to standardize a junction sequence.
        For more detail, please refer to the `source code <https://github.com/yutanagano/tidytcells>`_.

        .. code-block:: none

            IF input sequence contains non-amino acid symbols:
                set standardization status to failed

            IF a j symbol/locus and species are provided:
                attempt to resolve the correct conserved trailing amino acid (W / F / C)

            IF trimming is set to True and species is provided:
                attempt to trim the sequence based on conserved V gene patterns and trailing amin acid

            IF input sequence does not start with C and end with the correct conserved W / F / C:
                IF strict is set to True:
                    set standardization status to failed
                ELSE:
                    add C to the beginning and W / F / C to the end of the input sequence as required
                    set standardization status to successful
            ELSE:
                set standardization status to successful

            IF standardization status is set to successful:
                RETURN standardized sequence

            ELSE:
                IF on_fail is set to "reject":
                    RETURN None
                IF on_fail is set to "keep":
                    RETURN original sequence
    """
    seq = Parameter(seq, "seq").throw_error_if_not_of_type(str).value
    locus = (
        Parameter(locus, "locus")
        .throw_error_if_not_one_of("TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL", "IG", "TR", None)
        .value
    )
    j_symbol = (
        Parameter(j_symbol, "j_symbol")
        .throw_error_if_not_of_type(str, optional=True)
        .value
    )
    species = (
        Parameter(species, "species")
        .set_default("homosapiens")
        .throw_error_if_not_of_type(str)
        .value
    )
    strict = (
        Parameter(strict, "strict")
        .set_default(False)
        .throw_error_if_not_of_type(bool)
        .value
    )
    j_strict = (
        Parameter(j_strict, "j_strict")
        .set_default(False)
        .throw_error_if_not_of_type(bool)
        .value
    )
    trimming = (
        Parameter(trimming, "trimming")
        .set_default(False)
        .throw_error_if_not_of_type(bool)
        .value
    )
    on_fail = (
        Parameter(on_fail, "on_fail")
        .set_default("reject")
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

    original_input = seq

    seq = aa.standardize(seq=seq, on_fail="reject", log_failures=log_failures)

    not_valid_amino_acid_sequence = seq is None
    if not_valid_amino_acid_sequence:
        if on_fail == "reject":
            return None
        return original_input

    conserved_aa = "F"
    junction_matching_regex = re.compile(r"^C[A-Z]{4,}[FW]$")
    species = _utils.clean_and_lowercase(species)

    if j_symbol:
        conserved_aa = get_conserved_aa(j_symbol=j_symbol, locus=locus, species=species, log_failures=log_failures)

        if conserved_aa is not None:
            junction_matching_regex = re.compile(r"^C[A-Z]{4,}[" + conserved_aa + "]$")

        if conserved_aa is None:
            if j_strict:
                if on_fail == "reject":
                    return None
                return original_input
            else:
                logger.info(f"J symbol conserved amino acid could not be determined for {j_symbol}, using F as default.")
                conserved_aa = "F"

    # if trimming:
    #     # junction_matching_regex contains the correct ending aa pattern: [FW] (default)
    #     #   or a specific F/W/C if this was determined based on J gene/locus
    #
    seq = process_junction(seq, junction_matching_regex, conserved_aa, locus, species, trimming=trimming, check_motifs=True)

    if not junction_matching_regex.match(seq):
        if strict:
            if log_failures:
                logger.warning(
                    f"Failed to standardize {original_input}: not a valid junction sequence."
                )

            if on_fail == "reject":
                return None

            return original_input

    #     if not junction_matching_regex.match(seq):
    #         seq =  "C" + seq + conserved_aa
    #
    # return seq


def standardise(*args, **kwargs) -> Optional[str]:
    """
    Alias for :py:func:`tidytcells.junction.standardize`.
    """
    return standardize(*args, **kwargs)
