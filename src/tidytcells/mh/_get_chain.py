import logging
import re
from typing import Optional
from tidytcells._resources import VALID_HOMOSAPIENS_MH
from tidytcells._utils import Parameter


logger = logging.getLogger(__name__)


ALPHA_MATCHING_REGEX = re.compile(r"HLA-([ABCEFG]|D[PQR]A)")
BETA_MATCHING_REGEX = re.compile(r"HLA-D[PQR]B|B2M")


def get_chain(
    symbol: Optional[str] = None,
    log_failures: Optional[bool] = None,
    gene: Optional[str] = None,
    suppress_warnings: Optional[bool] = None,
) -> str:
    """
    Given a standardized MH gene / allele symbol, detect whether it codes for an alpha chain, beta chain, or beta-2 microglobulin (B2M) molecule.

    .. note::

        This function currently only recognises HLA (human leucocyte antigen or *Homo sapiens* MH), and not MH from other species.

    :param symbol:
        Standardized MH gene / allele symbol
    :type symbol:
        str
    :param log_failures:
        Report standardisation failures through logging (at level ``WARNING``).
        Defaults to ``True``.
    :type log_failures:
        bool
    :param gene:
        Alias for `symbol`.
    :type gene:
        str
    :param suppress_warnings:
        Disable warnings that are usually logged when standardisation fails.
        Deprecated in favour of `log_failures`.
    :type suppress_warnings:
        bool

    :return:
        ``'alpha'`` or ``'beta'`` if `symbol` is recognised and its chain is known, else ``None``.
    :rtype:
        Union[str, None]

    .. topic:: Example usage

        >>> tt.mh.get_chain("HLA-A")
        'alpha'
        >>> tt.mh.get_chain("HLA-DRB2")
        'beta'
        >>> tt.mh.get_chain("B2M")
        'beta'
    """
    symbol = (
        Parameter(symbol, "symbol")
        .resolve_with_alias(gene, "gene")
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

    symbol = symbol.split("*")[0]

    if not symbol in (*VALID_HOMOSAPIENS_MH, "B2M"):
        if log_failures:
            logger.warning(f"Unrecognized gene {symbol}. Is this standardized?")
        return None

    if ALPHA_MATCHING_REGEX.match(symbol):
        return "alpha"

    if BETA_MATCHING_REGEX.match(symbol):
        return "beta"

    if log_failures:
        logger.warning(f"Chain for {symbol} unknown.")

    return None
