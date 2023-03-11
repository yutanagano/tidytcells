import re
from .._utils.abstract_functions import standardise_aa_template
from warnings import warn


def standardise(seq: str, strict: bool = False, suppress_warnings: bool = False):
    """
    Ensures that a string value looks like a valid junction (CDR3) amino acid sequence.
    This function is a special variant of :py:func:`tidytcells.aa.standardise`.

    A valid junction sequence must:

    1. Be a valid amino acid sequence
    2. Begin with a cysteine (C)
    3. End with a phenylalanine (F) or a tryptophan (W)

    :param seq:
        String value representing a junction sequence.
    :type seq:
        ``str``
    :param strict:
        If ``True``, any string that does not look like a junction sequence is rejected.
        If ``False``, any inputs that are valid amino acid sequences but do not start with C and end with F/W are not rejected and instead are corrected by having a C appended to the beginning and an F appended at the end.
        Defaults to ``False``.
    :type strict:
        ``bool``
    :param suppress_warnings:
        Disable warnings that are usually emitted when standardisation fails.
        Defaults to ``False``.
    :type suppress_warnings:
        ``bool``

    :return:
        If possible, a standardised version of the input string is returned.
        If the input string cannot be standardised, it is rejected and ``None`` is returned.
    :rtype:
        ``str`` or ``None``

    .. topic:: Example usage

        Strings that look like junction sequences will be accepted, and returned in capitalised form.

        >>> tt.junction.standardise("csadaff")
        'CSADAFF'

        Strings that are valid amino acid sequences but do not stard and end with the appropriate residues will have a C and an F appended to its beginning and end respectively.

        >>> tt.junction.standardise("sadaf")
        'CSADAFF'

        However, setting ``strict`` to ``True`` will cause these cases to be rejected.

        >>> result = tt.junction.standardise("sadaf", strict=True)
        UserWarning: Input sadaf was rejected as it is not a valid junction sequence.
        >>> print(result)
        None
    """

    # take note of original input
    original_input = seq

    seq = standardise_aa_template(seq, suppress_warnings)

    if seq is None:  # not a valid amino acid sequence
        return None

    if not re.match(f"^C[A-Z]*[FW]$", seq):
        if strict:
            if not suppress_warnings:
                warn(
                    f"Input {original_input} was rejected as it is not a valid junction sequence."
                )
            return None
        seq = "C" + seq + "F"

    return seq
