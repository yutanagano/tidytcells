import re
from typing import Optional
import warnings

from tidytcells._resources import VALID_HOMOSAPIENS_MH
from tidytcells._utils import Parameter


ALPHA_MATCHING_REGEX = re.compile(r"HLA-([ABCEFG]|D[PQR]A)")
BETA_MATCHING_REGEX = re.compile(r"HLA-D[PQR]B|B2M")


def get_chain(
    gene: Optional[str] = None,
    suppress_warnings: bool = False,
) -> str:
    """
    Given a standardized MH gene name, detect whether it codes for an alpha chain, beta chain, or beta-2 microglobulin (B2M) molecule.

    .. note::

        This function currently only recognises HLA (human leucocyte antigen or *Homo sapiens* MH), and not MH from other species.

    :param gene:
        Standardized MH gene name
    :type gene:
        str
    :param suppress_warnings:
        Disable warnings that are usually emitted when chain classification fails.
        Defaults to ``False``.
    :type suppress_warnings:
        bool

    :return:
        ``'alpha'`` or ``'beta'`` if ``gene`` is recognised and its chain is known, else ``None``.
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
    Parameter(gene, "gene").throw_error_if_not_of_type(str)

    gene = gene.split("*")[0]

    if not gene in (*VALID_HOMOSAPIENS_MH, "B2M"):
        if not suppress_warnings:
            warnings.warn(f"Unrecognised gene {gene}. Is this standardized?")
        return None

    if ALPHA_MATCHING_REGEX.match(gene):
        return "alpha"

    if BETA_MATCHING_REGEX.match(gene):
        return "beta"

    if not suppress_warnings:
        warnings.warn(f"Chain for {gene} unknown.")
    return None
