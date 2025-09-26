import logging
from tidytcells._resources import (
    HOMOSAPIENS_TR_AA_SEQUENCES,
    HOMOSAPIENS_IG_AA_SEQUENCES,
    MUSMUSCULUS_TR_AA_SEQUENCES,
)
from typing import Optional, Tuple

TR_AA_SEQUENCES = {
    "homosapiens": HOMOSAPIENS_TR_AA_SEQUENCES,
    "musmusculus": MUSMUSCULUS_TR_AA_SEQUENCES,
}
IG_AA_SEQUENCES = {"homosapiens": HOMOSAPIENS_IG_AA_SEQUENCES}


logger = logging.getLogger(__name__)


def get_conserved_aa_for_j_symbol_for_species(
    j_symbol, species, log_failures
) -> Tuple[str, bool]:
    """
    Given a standardized J symbol and species, attempt to infer what the
    conserved residue at position 118 is.

    Returns:
        aa_118_target (str): The target (expected) residue at position 118.
            This is either the conserved residue inferred from the J symbol, or
            F if the conserved residue could not be inferred.
        aa_118_certain (bool): True if the conserved residue could be inferred
            with certainty, false otherwise.

    Raises:
        ValueError: If `j_symbol` is not recognized.
    """

    aa_dict = None

    if j_symbol.startswith("TR"):
        if species not in TR_AA_SEQUENCES:
            raise ValueError(
                f"Failed to determine the conserved residue at position 118 for J symbol {j_symbol}. "
                f"This feature is not supported for TR genes of species {species}."
            )
        aa_dict = TR_AA_SEQUENCES[species]
    elif j_symbol.startswith("IG"):
        if species not in IG_AA_SEQUENCES:
            raise ValueError(
                f"Failed to determine the conserved residue at position 118 for J symbol {j_symbol}. "
                f"This feature is not supported for IG genes of species {species}."
            )
        aa_dict = IG_AA_SEQUENCES[species]
    else:
        raise ValueError(
            f"Unrecognized J symbol {j_symbol}. "
            "Have you used tt.tr.standardize or tt.ig.standardize to standardize the symbol?"
        )

    if j_symbol in aa_dict:
        conserved_aa = _get_aa118_for_allele(aa_dict[j_symbol])
        if conserved_aa is not None:
            return (conserved_aa, True)
        return ("F", False)

    compatible_symbols = _get_compatible_symbols(j_symbol, aa_dict)
    possible_conserved_aas = set()

    for j_symbol_extended in compatible_symbols:
        conserved_aa = _get_aa118_for_allele(aa_dict[j_symbol_extended])
        if conserved_aa is not None:
            possible_conserved_aas.add(conserved_aa)

    if len(possible_conserved_aas) == 1:
        return (possible_conserved_aas.pop(), True)

    if len(possible_conserved_aas) == 0:
        if log_failures:
            logger.warning(
                f"No information available for conserved residue at position 118 for J symbol {j_symbol}. "
                f"Compatible symbols found: ({compatible_symbols})."
            )
        return ("F", False)

    if log_failures:
        logger.warning(
            f"The conserved residue at position 118 was ambiguous for J symbol {j_symbol}. "
            f"Compatible symbols found: ({compatible_symbols}). "
            f"Possible conserved residues: ({possible_conserved_aas})."
        )
    return ("F", False)


def _get_aa118_for_allele(allele_entry) -> Optional[str]:
    if "J-PHE" in allele_entry:
        return "F"

    if "J-TRP" in allele_entry:
        return "W"

    if "J-CYS" in allele_entry:
        return "C"

    return None


def _get_compatible_symbols(j_symbol, aa_dict):
    return [
        candidate
        for candidate in aa_dict.keys()
        if "J" in candidate
        and candidate.startswith(j_symbol)
        and not (j_symbol[-1].isnumeric() and candidate[len(j_symbol)].isnumeric())
    ]
