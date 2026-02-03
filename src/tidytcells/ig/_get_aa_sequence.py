from tidytcells._resources import SUPPORTED_SPECIES_AND_THEIR_IG_AA_SEQUENCES
from tidytcells import _utils
from tidytcells._utils import Parameter
from typing import Dict, Optional




def get_aa_sequence(
    symbol: Optional[str] = None,
    species: Optional[str] = None,
    gene: Optional[str] = None,
) -> Dict[str, str]:
    """
    Look up the amino acid sequence of a given IG allele.

    .. topic:: Supported species

        - ``"homosapiens"``

    :param symbol:
        Standardized allele symbol.
        Note that the symbol must be specified to the level of the allele.
        Note that some alleles, notably those of non-functional genes, will not have resolvable amino acid sequences.
    :type symbol:
        str
    :param species:
        Species to which the IG gene in question belongs (see above for supported species).
        Defaults to ``"homosapiens"``.
    :type species:
        str
    :param gene:
        Alias for `symbol`.
    :type gene:
        str

    :return:
        A dictionary with keys corresponding to names of different sequence regions within the allele, and values corresponding to their amino acid sequences.
    :rtype:
        Dict[str, str]

    .. topic:: Example usage

        Get amino acid sequence information about the human V gene IGHV1-18*02.

        >>> tt.ig.get_aa_sequence(gene="IGHV1-18*02", species="homosapiens")
        {'CDR1-IMGT': 'GYTFTSYG', 'CDR2-IMGT': 'ISAYNGNT', 'FR1-IMGT': 'QVQLVQSGAEVKKPGASVKVSCKAS', 'FR2-IMGT': 'ISWVRQAPGQGLEWMGW', 'FR3-IMGT': 'NYAQKLQGRVTMTTDTSTSTAYMELRSLRSDDTA', 'V-REGION': 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYGISWVRQAPGQGLEWMGWISAYNGNTNYAQKLQGRVTMTTDTSTSTAYMELRSLRSDDTA'}

        Get amino acid sequence information about the human J gene IGLJ1*01.

        >>> tt.ig.get_aa_sequence(gene="IGLJ1*01", species="homosapiens")
        {'J-REGION': 'YVFGTGTKVTVL'}

    """
    symbol = (
        Parameter(symbol, "symbol")
        .resolve_with_alias(gene, "gene")
        .throw_error_if_not_of_type(str)
        .value
    )
    species = (
        Parameter(species, "species")
        .set_default("homosapiens")
        .throw_error_if_not_of_type(str)
        .value
    )

    species = _utils.clean_and_lowercase(species)

    species_is_supported = species in SUPPORTED_SPECIES_AND_THEIR_IG_AA_SEQUENCES
    if not species_is_supported:
        raise ValueError(f"Unsupported species: {species}. No data available.")

    aa_sequence_dict = SUPPORTED_SPECIES_AND_THEIR_IG_AA_SEQUENCES[species]

    if symbol in aa_sequence_dict:
        return aa_sequence_dict[symbol]

    raise ValueError(f"No data found for IG gene {symbol} for species {species}.")
