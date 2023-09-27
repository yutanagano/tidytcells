import re
from typing import Optional
import warnings

from tidytcells._resources import VALID_HOMOSAPIENS_MH
from tidytcells._utils import Parameter


CLASS_1_MATCHING_REGEX = re.compile(r"HLA-[ABCEFG]|B2M")
CLASS_2_MATCHING_REGEX = re.compile(r"HLA-D[PQR][AB]")


def get_class(
    gene: Optional[str] = None,
    suppress_warnings: bool = False,
) -> int:
    """
    Given a standardized MH gene name, detect whether it comprises a class I (MH1) or II (MH2) receptor.

    .. note::

        This function currently only recognises HLA (human leucocyte antigen or *Homo sapiens* MH), and not MH from other species.

    :param gene:
        Standardized MH gene name
    :type gene:
        str
    :param suppress_warnings:
        Disable warnings that are usually emitted when classification fails.
        Defaults to ``False``.
    :type suppress_warnings:
        bool

    :return:
        ``1`` or ``2`` if ``gene`` is recognised and its class is known, else ``None``.
    :rtype:
        Union[int, None]

    .. topic:: Example usage

        >>> tt.mh.get_class("HLA-A")
        1
        >>> tt.mh.get_class("HLA-DRB2")
        2
        >>> tt.mh.get_class("B2M")
        1
    """
    Parameter(gene, "gene").throw_error_if_not_of_type(str)

    gene = gene.split("*")[0]

    if not gene in (*VALID_HOMOSAPIENS_MH, "B2M"):
        if not suppress_warnings:
            warnings.warn(f"Unrecognised gene {gene}. Is this standardized?")
        return None

    if CLASS_1_MATCHING_REGEX.match(gene):
        return 1

    if CLASS_2_MATCHING_REGEX.match(gene):
        return 2

    if not suppress_warnings:
        warnings.warn(f"Class for {gene} unknown.")
    return None
