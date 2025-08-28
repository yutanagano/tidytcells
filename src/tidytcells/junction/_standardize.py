import logging
import re
from tidytcells import aa, _utils
from typing import Literal, Optional
from tidytcells._utils.parameter import Parameter
from tidytcells._utils.conserved_aa_lookup import get_conserved_aa_for_j_symbol_for_species

logger = logging.getLogger(__name__)



def standardize(
    seq: str,
    j_symbol: Optional[str] = None,
    species: Optional[str] = None,
    strict: Optional[bool] = None,
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
    3. End with a phenylalanine (F) or a tryptophan (W)

    :param seq:
        String value representing a junction sequence.
    :type seq:
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

            IF a j symbol and species are provided:
                attempt to resolve the correct conserved trailing amino acid (W / F)

            IF input sequence does not start with C and end with the correct conserved W / F:
                IF strict is set to True:
                    set standardization status to failed
                ELSE:
                    add C to the beginning and W / F to the end of the input sequence as required
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
    junction_matching_regex = re.compile(r"^C[A-Z]*[FW]$")

    if j_symbol:
        species = _utils.clean_and_lowercase(species)
        conserved_aa = get_conserved_aa_for_j_symbol_for_species(j_symbol, species, log_failures=log_failures)
        if conserved_aa is None:
            if on_fail == "reject":
                return None
            return original_input
        
        junction_matching_regex = re.compile(fr"^C[A-Z]*{conserved_aa}$")

    if not junction_matching_regex.match(seq):
        if strict:
            if log_failures:
                logger.warning(
                    f"Failed to standardize {original_input}: not a valid junction sequence."
                )

            if on_fail == "reject":
                return None

            return original_input

        if not seq.startswith("C"):
            seq = "C" + seq

        if not junction_matching_regex.match(seq):
            seq = seq + conserved_aa

    return seq


def standardise(*args, **kwargs) -> Optional[str]:
    """
    Alias for :py:func:`tidytcells.junction.standardize`.
    """
    return standardize(*args, **kwargs)
