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
    gene_name: Optional[str] = None,
) -> str:
    """
    Given a standardized MH gene name, detect whether it codes for an alpha or a beta chain molecule.

    .. note::

        This function currently only recognises HLA, and not MH from other species.

    :param gene:
        Standardized MH gene name
    :type gene:
        str
    :param suppress_warnings:
        Disable warnings that are usually emitted when chain classification fails.
        Defaults to ``False``.
    :type suppress_warnings:
        bool

    :param gene_name:
        Alias for the parameter ``gene``.

        .. caution:: This will be deprecated soon in favour of ``gene``.
    :type gene_name:
        str

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
    gene = Parameter(gene, "gene").resolve_with_alias_and_return_value(
        Parameter(gene_name, "gene_name")
    )
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
