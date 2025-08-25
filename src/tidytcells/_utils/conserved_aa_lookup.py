import logging

from tidytcells._resources import HOMOSAPIENS_TR_AA_SEQUENCES, HOMOSAPIENS_IG_AA_SEQUENCES, MUSMUSCULUS_TR_AA_SEQUENCES


logger = logging.getLogger(__name__)



def _get_conserved_aa_exact_symbol(aa_dict, j_symbol):
    if "J-PHE" in aa_dict[j_symbol] and aa_dict[j_symbol]["J-PHE"] == "F":
        return "F"
    elif "J-TRP" in aa_dict[j_symbol] and aa_dict[j_symbol]["J-TRP"] == "W":
        return "W"


def _is_valid_extension(original, extension):
    if extension.startswith(original):
        if original[-1].isnumeric() and extension[len(original)].isnumeric():
            return False # cannot concatenate numbers if original ends with a number (e.g., TRAV1 to TRAV13 is not valid)
        return True
    return False


def _get_all_extended_symbols(j_symbol, aa_dict):
    return sorted(list({key for key in aa_dict.keys() if _is_valid_extension(j_symbol, key)}))


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


def _get_conserved_aa(j_symbol, aa_dict, log_failures):
    if j_symbol in aa_dict:
        return _get_conserved_aa_exact_symbol(aa_dict, j_symbol)

    return _resolve_conserved_aa_from_partial_j_symbol(j_symbol, aa_dict, log_failures)


def _get_aa_dict(gene_type, species):
    if gene_type == "TR" and species == "homosapiens":
        return HOMOSAPIENS_TR_AA_SEQUENCES
    elif gene_type == "IG" and species == "homosapiens":
        return HOMOSAPIENS_IG_AA_SEQUENCES
    elif gene_type == "TR" and species == "musmusculus":
        return MUSMUSCULUS_TR_AA_SEQUENCES


def get_conserved_aa_for_j_symbol_for_species(j_symbol, species, log_failures):
    if j_symbol.startswith("TR") or j_symbol.startswith("IG"):
        aa_dict = _get_aa_dict(j_symbol[0:2], species)

        if aa_dict is not None:
            return _get_conserved_aa(j_symbol, aa_dict, log_failures)

        else:
            raise ValueError(f"Failed to determine the conserved trailing amino acid for J symbol {j_symbol}: this feature "
                f"is not supported for {j_symbol[0:2]} genes for species {species}.")
    else:
        raise ValueError(f"Failed to determine the conserved trailing amino acid for J symbol {j_symbol}: symbol is not formatted correctly."
                         f"Please use tt.tr.standardize or tt.ig.standardize to correct the symbol.")
