import logging
import re
from tidytcells import aa, _utils
from typing import Literal, Optional
from tidytcells._utils.parameter import Parameter
from tidytcells._utils.conserved_aa_lookup import (
    get_conserved_aa_for_j_symbol_for_species,
)

logger = logging.getLogger(__name__)


def standardize(
    seq: str,
    j_symbol: Optional[str] = None,
    species: Optional[str] = None,
    allow_uncertain_118: Optional[bool] = None,
    add_missing_conserved: Optional[bool] = None,
    on_fail: Optional[Literal["reject", "keep"]] = None,
    log_failures: Optional[bool] = None,
    j_strict: Optional[bool] = None,
    strict: Optional[bool] = None,
    suppress_warnings: Optional[bool] = None,
) -> Optional[str]:
    """
    Ensures that a string value looks like a valid junction (CDR3) amino acid
    sequence.

    A valid junction sequence must:

    1. Be a valid amino acid sequence
    2. Begin with a cysteine (C)
    3. End with a phenylalanine (F), tryptophan (W) or cysteine (C) in a way
       consistent with `j_symbol` if supplied

    :param seq:
        The junction sequence.
    :type seq:
        str
    :param j_symbol:
        The TR/IG J symbol used to determine the correct conserved trailing
        amino acid at position 118 (F / W / C). If the symbol does not resolve
        to a single allele but all productive alleles consistent with the
        symbol have the same conserved residue, this will be set as the
        expected ending residue. If the supplied symbol does not map to any
        (group of) known J alleles, the function will raise a ``ValueError``.
    :type j_symbol:
        str
    :param species:
        The species that produced the underlying receptor. Defaults to
        ``homosapiens``.
    :type species:
        str
    :param allow_uncertain_118:
        If ``False``, standardization immediately fails if the expected
        conserved trailing amino acid at position 118 cannot be determined with
        certainty using `j_symbol`, or if `j_symbol` is not supplied. If
        ``True``, in the event of an uncertain residue at position 118, either
        F or W is accepted, and if a trailing residue must be appended (see
        parameter `add_missing_conserved`), an F will be added. Defaults to
        ``True``.
    :type allow_uncertain_118:
        bool
    :param add_missing_conserved:
        If ``False``, standardization immediately fails for any input sequence
        that does not start and end with the expected conserved residues. If
        ``True``, any inputs that are valid amino acid sequences but do not
        start and end as expected are corrected by adding a C at the beginning
        and the expected trailing residue (see `allow_uncertain_118`) at the
        end. Defaults to ``True``.
    :type add_missing_conserved:
        bool
    :param on_fail:
        Behaviour when standardization fails. If set to ``"reject"``, returns
        ``None`` on failure. If set to ``"keep"``, returns the original input.
        Defaults to ``"reject"``.
    :type on_fail:
        str
    :param log_failures:
        Report standardisation failures through logging (at level ``WARNING``).
        Defaults to ``True``.
    :type log_failures:
        bool
    :param j_strict:
        Inverse setting to `allow_uncertain_118`. Deprecated in favor of
        `allow_uncertain_118`.
    :type j_strict:
        bool
    :param strict:
        Inverse setting to `add_missing_conserved`. Deprecated in favor of
        `add_missing_conserved`.
    :type strict:
        bool
    :param suppress_warnings:
        Disable warnings that are usually logged when standardisation fails.
        Deprecated in favour of `log_failures`.
    :type suppress_warnings:
        bool

    :return:
        If possible, a standardized version of the input string is returned. If
        the input string cannot be standardized, the function follows the
        behaviour as set by `on_fail`.
    :rtype:
        Optional[str]

    .. topic:: Example usage

        Strings that look like junction sequences will be accepted, and
        returned in capitalised form.

        >>> tt.junction.standardize("csadaf")
        'CSADAF'

        Strings that are valid amino acid sequences but do not start and end
        with the appropriate residues will have a C and the appropriate
        conserved trailing residue at position 118 (defaults to F) appended to
        its beginning and end as required.

        >>> tt.junction.standardize("sada")
        'CSADAF'

        The conserved trailing residue can be intelligently inferred if
        `j_symbol` is supplied.

        >>> tt.junction.standardize("sada", j_symbol="TRAJ38*01")
        'CSADAW'

        Furthermore, setting `add_missing_conserved` to ``False`` will cause
        these cases to be rejected.

        >>> result = tt.junction.standardize("sada", add_missing_conserved=False)
        Input sadaf was rejected as it is not a valid junction sequence.
        >>> print(result)
        None

    .. topic:: Decision Logic

        To provide an easy way to gauge the scope and limitations of
        standardization, below is a simplified overview of the decision logic
        employed when attempting to standardize a junction sequence. For more
        detail, please refer to the
        `source code <https://github.com/yutanagano/tidytcells>`_.

        .. code-block:: none

            IF input sequence contains non-amino acid symbols:
                set standardization status to failed
                skip rest of standardization

            // inferred using J symbol if supplied
            IF expected trailing residue at position 118 uncertain:
                {
                    IF allow_uncertain_118:
                        accept either F or W
                    ELSE:
                        set standardization status to failed
                        skip rest of standardization
                }

            IF input sequence starts (C) and ends (F / W / C) as expected:
                set standardization status to successful
            ELSE:
                {
                    IF add_missing_conserved:
                        append expected starting and ending residues
                        set standardization status to successful
                    ELSE:
                        set standardization status to failed
                }

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
    j_strict_inverted = not j_strict if j_strict is not None else None
    allow_uncertain_118 = (
        Parameter(allow_uncertain_118, "allow_uncertain_118")
        .set_default(True)
        .resolve_with_alias(j_strict_inverted, "j_strict")
        .throw_error_if_not_of_type(bool)
        .value
    )
    strict_inverted = not strict if strict is not None else None
    add_missing_conserved = (
        Parameter(add_missing_conserved, "add_missing_conserved")
        .set_default(True)
        .resolve_with_alias(strict_inverted, "strict")
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

    if seq is None:
        if on_fail == "reject":
            return None
        return original_input

    aa_118_target = "F"
    aa_118_certain = False
    junction_matching_regex = None

    if j_symbol:
        species = _utils.clean_and_lowercase(species)
        aa_118_target, aa_118_certain = get_conserved_aa_for_j_symbol_for_species(
            j_symbol, species, log_failures=log_failures
        )

    if aa_118_certain:
        junction_matching_regex = re.compile(rf"^C[A-Z]*{aa_118_target}$")
    else:
        if not allow_uncertain_118:
            if on_fail == "reject":
                return None
            return original_input
        else:
            logger.info(
                f"Unclear residue at position 118 (j_symbol = {j_symbol}), accepting either F or W."
            )
            junction_matching_regex = re.compile(r"^C[A-Z]*[FW]$")

    if not junction_matching_regex.match(seq):
        if not add_missing_conserved:
            if log_failures:
                logger.warning(
                    f"Failed to standardize {original_input}: not a valid junction sequence."
                )

            if on_fail == "reject":
                return None

            return original_input

        seq = "C" + seq + aa_118_target

    return seq


def standardise(*args, **kwargs) -> Optional[str]:
    """
    Alias for :py:func:`tidytcells.junction.standardize`.
    """
    return standardize(*args, **kwargs)
