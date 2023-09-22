import warnings

from tidytcells._resources import AMINO_ACIDS
from tidytcells._utils import Parameter


def standardize(seq: str, on_fail: str = "reject", suppress_warnings: bool = False):
    """
    Ensures that a string value looks like a valid amino acid sequence.

    :param seq:
        String value representing an amino acid sequence.
    :type seq:
        str
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
        Capitalised version of ``seq``, if seq is a valid amino acid sequence.
        Otherwise follow behaviour set by ``on_fail``.
    :rtype:
        Union[str, None]

    .. topic:: Example usage

        Strings that look like amino acid sequences will be accepted, and returned in capitalised form.

        >>> tt.aa.standardize("sqllnakyl")
        'SQLLNAKYL'

        Any strings that contain characters that cannot be recognised as amino acids will be rejected, and the function will return ``None``.

        >>> result = tt.aa.standardize("sqll?akyl")
        UserWarning: Input sqll?akyl was rejected as it is not a valid amino acid sequence.
        >>> print(result)
        None

    .. topic:: Decision Logic

        To provide an easy way to gauge the scope and limitations of standardization, below is a simplified overview of the decision logic employed when attempting to standardize an amino acid sequence.
        For more detail, please refer to the `source code <https://github.com/yutanagano/tidytcells>`_.

        .. code-block:: none

            IF input sequence contains non-amino acid symbols:
                set standardization status to failed
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
    Parameter(seq, "seq").throw_error_if_not_of_type(str)
    Parameter(on_fail, "on_fail").throw_error_if_not_one_of("reject", "keep")
    Parameter(suppress_warnings, "suppress_warnings").throw_error_if_not_of_type(bool)

    original_input = seq

    seq = seq.upper()

    for char in seq:
        if not char in AMINO_ACIDS:
            if not suppress_warnings:
                warnings.warn(
                    f"Failed to standardize {original_input}: not a valid amino acid sequence."
                )
            if on_fail == "reject":
                return None
            return original_input

    return seq


def standardise(*args, **kwargs):
    """
    Alias for :py:func:`tidytcells.aa.standardize`.
    """
    return standardize(*args, **kwargs)
