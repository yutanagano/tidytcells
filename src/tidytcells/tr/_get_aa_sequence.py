from typing import Dict

from tidytcells._resources import (
    HOMOSAPIENS_TR_AA_SEQUENCES,
    MUSMUSCULUS_TR_AA_SEQUENCES,
)
from tidytcells import _utils
from tidytcells._utils import Parameter


SUPPORTED_SPECIES_AND_THEIR_AA_SEQUENCES = {
    "homosapiens": HOMOSAPIENS_TR_AA_SEQUENCES,
    "musmusculus": MUSMUSCULUS_TR_AA_SEQUENCES,
}


def get_aa_sequence(gene: str, species: str = "homosapiens") -> Dict[str, str]:
    """
    Look up the amino acid sequence of a given TR gene.

    .. topic:: Supported species

        - ``"homosapiens"``
        - ``"musmusculus"``

    :param gene:
        Standardized gene name.
        The gene must be specified to the level of the allele.
        Note that some genes, notably the non-functional ones, will not have resolvable amino acid sequences.
    :type gene:
        str
    :param species:
        Species to which the TR gene in question belongs (see above for supported species).
        Defaults to ``"homosapiens"``.
    :type species:
        str

    :return:
        A dictionary with keys corresponding to names of different sequence regions within the gene, and values corresponding to their amino acid sequences.
    :rtype:
        Dict[str, str]

    .. topic:: Example usage

        Get amino acid sequence information about the human V gene TRBV2*01.

        >>> tt.tr.get_aa_sequence(gene="TRBV2*01", species="homosapiens")
        {'CDR1-IMGT': 'SNHLY', 'CDR2-IMGT': 'FYNNEI', 'FR1-IMGT': 'EPEVTQTPSHQVTQMGQEVILRCVPI', 'FR2-IMGT': 'FYWYRQILGQKVEFLVS', 'FR3-IMGT': 'SEKSEIFDDQFSVERPDGSNFTLKIRSTKLEDSAMYFC', 'V-REGION': 'EPEVTQTPSHQVTQMGQEVILRCVPISNHLYFYWYRQILGQKVEFLVSFYNNEISEKSEIFDDQFSVERPDGSNFTLKIRSTKLEDSAMYFCASSE'}

        Get amino acid sequence information about the *Mus musculus* J gene TRAJ32*02.

        >>> tt.tr.get_aa_sequence(gene="TRAJ32*02", species="musmusculus")
        {'FR4-IMGT': 'FGTGTLLSVKP', 'J-REGION': 'NYGGSGNKLIFGTGTLLSVKP'}
    """
    Parameter(gene, "gene").throw_error_if_not_of_type(str)
    Parameter(species, "species").throw_error_if_not_of_type(str)

    species = _utils.clean_and_lowercase(species)

    species_is_supported = species in SUPPORTED_SPECIES_AND_THEIR_AA_SEQUENCES
    if not species_is_supported:
        raise ValueError(f"Unsupported species: {species}. No data available.")

    aa_sequence_dict = SUPPORTED_SPECIES_AND_THEIR_AA_SEQUENCES[species]

    if gene in aa_sequence_dict:
        return aa_sequence_dict[gene]

    raise ValueError(f"No data found for TR gene {gene} for species {species}.")
