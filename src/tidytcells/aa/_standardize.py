import logging
from tidytcells._resources import AMINO_ACIDS
from tidytcells._utils import Parameter
from typing import Optional, Literal


logger = logging.getLogger(__name__)


def standardize(
    seq: str,
    on_fail: Optional[Literal["reject", "keep"]] = None,
    log_failures: Optional[bool] = None,
    suppress_warnings: Optional[bool] = None,
):
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
        Capitalised version of `seq`, if seq is a valid amino acid sequence.
        Otherwise follow behaviour set by `on_fail`.
    :rtype:
        Union[str, None]

    .. topic:: Example usage

        Strings that look like amino acid sequences will be accepted, and returned in capitalised form.

        >>> tt.aa.standardize("sqllnakyl")
        'SQLLNAKYL'

        Any strings that contain characters that cannot be recognised as amino acids will be rejected, and the function will return ``None``.

        >>> result = tt.aa.standardize("sqll?akyl")
        Input sqll?akyl was rejected as it is not a valid amino acid sequence.
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
    seq = Parameter(seq, "seq").throw_error_if_not_of_type(str).value
    on_fail = (
        Parameter(on_fail, "on_fail")
        .set_default("reject")
        .throw_error_if_not_one_of("reject", "keep")
        .value
    )
    log_failures = (
        Parameter(log_failures, "suppress_warnings")
        .set_default(True)
        .resolve_with_alias(suppress_warnings, "suppress_warnings")
        .throw_error_if_not_of_type(bool)
        .value
    )

    original_input = seq

    seq = seq.upper()

    for char in seq:
        if not char in AMINO_ACIDS:
            if log_failures:
                logger.warning(
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
