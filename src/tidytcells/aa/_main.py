from .._utils.abstract_functions import standardise_aa_template


def standardise(seq: str, suppress_warnings: bool = False):
    """
    Ensures that a string value looks like a valid amino acid sequence.

    :param seq:
        String value representing an amino acid sequence.
    :type seq:
        ``str``
    :param suppress_warnings:
        Disable warnings that are usually emitted when standardisation fails.
        Defaults to ``False``.
    :type suppress_warnings:
        ``bool``

    :return:
        Capitalised version of ``seq``, if seq is a valid amino acid sequence.
        Otherwise the input is rejected and ``None`` is returned.
    :rtype:
        ``str`` or ``None``

    .. topic:: Example usage

        Strings that look like amino acid sequences will be accepted, and returned in capitalised form.

        >>> tt.aa.standardise("sqllnakyl")
        'SQLLNAKYL'

        Any strings that contain characters that cannot be recognised as amino acids will be rejected, and the function will return ``None``.

        >>> result = tt.aa.standardise("sqll?akyl")
        UserWarning: Input sqll?akyl was rejected as it is not a valid amino acid sequence.
        >>> print(result)
        None
    """

    return standardise_aa_template(seq, suppress_warnings)
