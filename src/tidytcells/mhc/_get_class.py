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
    gene_name: Optional[str] = None,
) -> int:
    """
    Given a standardized MHC gene name, detect whether it comprises a class I
    or II MHC receptor complex.

    .. note::

        This function currently only recognises HLAs, and not MHCs from other species.

    :param gene:
        Standardized MHC gene name
    :type gene:
        str
    :param suppress_warnings:
        Disable warnings that are usually emitted when classification fails.
        Defaults to ``False``.
    :type suppress_warnings:
        bool

    :param gene_name:
        Alias for the parameter ``gene``.

        .. caution:: This will be deprecated soon in favour of ``gene``.
    :type gene_name:
        str

    :return:
        ``1`` or ``2`` if ``gene`` is recognised and its class is known, else ``None``.
    :rtype:
        Union[int, None]

    .. topic:: Example usage

        >>> tt.mhc.get_class("HLA-A")
        1
        >>> tt.mhc.get_class("HLA-DRB2")
        2
        >>> tt.mhc.get_class("B2M")
        1
    """
    gene = Parameter(gene, "gene").resolve_with_alias_and_return_value(
        Parameter(gene_name, "gene_name")
    )
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


def classify(*args, **kwargs):
    """
    .. caution:: This will be deprecated soon in favour of :py:func:`tidytcells.mhc.get_class`.

    Alias for :py:func:`tidytcells.mhc.get_class`.
    """
    warnings.warn(
        '"mhc.classify" as an alias will be deprecated in the near future. Please switch to using "mhc.get_class".',
        FutureWarning,
    )
    return get_class(*args, **kwargs)
