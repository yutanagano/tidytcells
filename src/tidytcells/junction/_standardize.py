import re
import warnings

from tidytcells import aa


JUNCTION_MATCHING_REGEX = re.compile(f"^C[A-Z]*[FW]$")


def standardize(
    seq: str,
    strict: bool = False,
    on_fail: str = "reject",
    suppress_warnings: bool = False,
):
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
    :param suppress_warnings:
        Disable warnings that are usually emitted when standardisation fails.
        Defaults to ``False``.
    :type suppress_warnings:
        bool

    :return:
        If possible, a standardized version of the input string is returned.
        If the input string cannot be standardized, the function follows the behaviour as set by ``on_fail``.
    :rtype:
        Union[str, None]

    .. topic:: Example usage

        Strings that look like junction sequences will be accepted, and returned in capitalised form.

        >>> tt.junction.standardize("csadaff")
        'CSADAFF'

        Strings that are valid amino acid sequences but do not stard and end with the appropriate residues will have a C and an F appended to its beginning and end respectively.

        >>> tt.junction.standardize("sadaf")
        'CSADAFF'

        However, setting ``strict`` to ``True`` will cause these cases to be rejected.

        >>> result = tt.junction.standardize("sadaf", strict=True)
        UserWarning: Input sadaf was rejected as it is not a valid junction sequence.
        >>> print(result)
        None

    .. topic:: Decision Logic

        To provide an easy way to gauge the scope and limitations of standardization, below is a simplified overview of the decision logic employed when attempting to standardize a junction sequence.
        For more detail, please refer to the `source code <https://github.com/yutanagano/tidytcells>`_.

        .. code-block:: none

            IF input sequence contains non-amino acid symbols:
                set standardization status to failed

            IF input sequence does not start with C and end with F:
                IF strict is set to True:
                    set standardization status to failed
                ELSE:
                    add C to the beginning and F to the end of the input sequence
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
    original_input = seq

    seq = aa.standardize(seq=seq, on_fail="reject", suppress_warnings=suppress_warnings)

    not_valid_amino_acid_sequence = seq is None
    if not_valid_amino_acid_sequence:
        if on_fail == "reject":
            return None
        return original_input

    if not JUNCTION_MATCHING_REGEX.match(seq):
        if strict:
            if not suppress_warnings:
                warnings.warn(
                    f"Failed to standardize {original_input}: not a valid junction sequence."
                )
            if on_fail == "reject":
                return None
            return original_input
        seq = "C" + seq + "F"

    return seq


def standardise(*args, **kwargs):
    """
    Alias for :py:func:`tidytcells.junction.standardize`.
    """
    return standardize(*args, **kwargs)
