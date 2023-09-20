from typing import Dict

from tidytcells._resources import HOMOSAPIENS_TCR_AA_SEQUENCES
from tidytcells import _utils
from tidytcells._utils import Parameter


def get_aa_sequence(gene: str, species: str = "homosapiens") -> Dict[str, str]:
    """
    Look up the amino acid sequence of a given TCR gene.

    .. topic:: Supported species

        - ``'homosapiens'``

    .. note::

        This function currently only supports V genes.
        Support for J genes is planned for the future.

    :param gene:
        Standardized gene name.
        The gene must be specified to the level of the allele.
        Note that some genes, notably the non-functional ones, will not have resolvable amino acid sequences.
    :type gene:
        str
    :param species:
        Species to which the TCR gene in question belongs (see above for supported species).
        Defaults to ``'homosapiens'``.
    :type species:
        str

    :return:
        A dictionary with keys corresponding to names of different sequence regions within the gene, and values corresponding to their amino acid sequences.
    :rtype:
        Dict[str, str]

    .. topic:: Example usage

        Get amino acid sequence information about the human V gene TRBV2*01.

        >>> tt.tcr.get_aa_sequence(gene="TRBV2*01", species="homosapiens")
        {'CDR1-IMGT': 'SNHLY', 'CDR2-IMGT': 'FYNNEI', 'FR1-IMGT': 'EPEVTQTPSHQVTQMGQEVILRCVPI', 'FR2-IMGT': 'FYWYRQILGQKVEFLVS', 'FR3-IMGT': 'SEKSEIFDDQFSVERPDGSNFTLKIRSTKLEDSAMYFC', 'V-REGION': 'EPEVTQTPSHQVTQMGQEVILRCVPISNHLYFYWYRQILGQKVEFLVSFYNNEISEKSEIFDDQFSVERPDGSNFTLKIRSTKLEDSAMYFCASSE'}
    """
    Parameter(gene, "gene").throw_error_if_not_of_type(str)
    Parameter(species, "species").throw_error_if_not_of_type(str)

    species = _utils.lowercase_and_remove_whitespace(species)

    # Currently only supports homosapiens
    if species != "homosapiens":
        raise ValueError(f"Unsupported species: {species}. No data available.")

    if gene in HOMOSAPIENS_TCR_AA_SEQUENCES:
        return HOMOSAPIENS_TCR_AA_SEQUENCES[gene]

    raise ValueError(f"No data found for TCR gene {gene} for species {species}.")
