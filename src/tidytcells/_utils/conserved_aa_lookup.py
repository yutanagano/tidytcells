import logging
from tidytcells._resources import (
    HOMOSAPIENS_TR_AA_SEQUENCES,
    HOMOSAPIENS_IG_AA_SEQUENCES,
    MUSMUSCULUS_TR_AA_SEQUENCES,
)

TR_AA_SEQUENCES = {
    "homosapiens": HOMOSAPIENS_TR_AA_SEQUENCES,
    "musmusculus": MUSMUSCULUS_TR_AA_SEQUENCES,
}
IG_AA_SEQUENCES = {"homosapiens": HOMOSAPIENS_IG_AA_SEQUENCES}


logger = logging.getLogger(__name__)


def get_conserved_aa_for_j_symbol_for_species(j_symbol, species, log_failures):
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

    if j_symbol.startswith("TR"):
        if species not in TR_AA_SEQUENCES:
            raise ValueError(
                f"Failed to determine the conserved trailing amino acid for J symbol {j_symbol}: this feature "
                f"is not supported for TR genes for species {species}."
            )
        return _get_conserved_aa(j_symbol, TR_AA_SEQUENCES[species], log_failures)

    if j_symbol.startswith("IG"):
        if species not in IG_AA_SEQUENCES:
            raise ValueError(
                f"Failed to determine the conserved trailing amino acid for J symbol {j_symbol}: this feature "
                f"is not supported for IG genes for species {species}."
            )
        return _get_conserved_aa(j_symbol, IG_AA_SEQUENCES[species], log_failures)

    raise ValueError(
        f"Unrecognized J symbol {j_symbol}. "
        "Have you used tt.tr.standardize or tt.ig.standardize to standardize the symbol?"
    )


def _get_conserved_aa(j_symbol, aa_dict, log_failures):
    if j_symbol in aa_dict:
        return _get_conserved_aa_exact_symbol(aa_dict, j_symbol)

    return _resolve_conserved_aa_from_partial_j_symbol(j_symbol, aa_dict, log_failures)


def _get_conserved_aa_exact_symbol(aa_dict, j_symbol):
    if "J-PHE" in aa_dict[j_symbol] and aa_dict[j_symbol]["J-PHE"] == "F":
        return "F"
    elif "J-TRP" in aa_dict[j_symbol] and aa_dict[j_symbol]["J-TRP"] == "W":
        return "W"
    elif "J-CYS" in aa_dict[j_symbol] and aa_dict[j_symbol]["J-CYS"] == "C":
        return "C"


def _is_valid_extension(original, extension):
    if extension.startswith(original):
        if original[-1].isnumeric() and extension[len(original)].isnumeric():
            return False  # cannot concatenate numbers if original ends with a number (e.g., TRAV1 to TRAV13 is not valid)
        return True
    return False


def _get_all_extended_symbols(j_symbol, aa_dict):
    return sorted(
        list(
            {
                key
                for key in aa_dict.keys()
                if _is_valid_extension(j_symbol, key)
                if "J" in key
            }
        )
    )


def _resolve_conserved_aa_from_partial_j_symbol(j_symbol, aa_dict, log_failures):
    extended_symbols = _get_all_extended_symbols(j_symbol, aa_dict)
    all_results = set()

    if len(extended_symbols) == 0:
        if log_failures:
            logger.warning(
                f"Failed to determine the conserved trailing amino acid for J symbol {j_symbol}: "
                f"symbol not recognized, and no extended allele names (starting with {j_symbol}) were found."
            )
        return None

    for j_symbol_extended in extended_symbols:
        new_aa = _get_conserved_aa_exact_symbol(aa_dict, j_symbol_extended)
        if new_aa is not None:
            all_results.add(new_aa)

    if len(all_results) == 1:
        return all_results.pop()
    elif len(all_results) > 1:
        if log_failures:
            logger.warning(
                f"Failed to determine the conserved trailing amino acid for J symbol {j_symbol}: "
                f"conserved amino acid was ambiguous based on extended allele names ({extended_symbols}). "
            )
        return None
