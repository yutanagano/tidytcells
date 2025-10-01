import logging
from util import remove_allele
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






def _get_aa118_for_allele(allele_entry) -> Optional[str]:
    if "J-PHE" in allele_entry:
        return "F"

    if "J-TRP" in allele_entry:
        return "W"

    if "J-CYS" in allele_entry:
        return "C"

    if "J-VAL" in allele_entry:
        return "V"

    return None


def _get_compatible_symbols(j_symbol, aa_dict, gene):
    return [
        candidate
        for candidate in aa_dict.keys()
        if gene in candidate
        and candidate.startswith(j_symbol)
        and not (j_symbol[-1].isnumeric() and candidate[len(j_symbol)].isnumeric())
           and not (j_symbol[-1].isnumeric() and candidate[len(j_symbol)] == "P")  # do not add 'pseudogene'
    ]



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

    compatible_symbols = _get_compatible_symbols(j_symbol, aa_dict, gene="J")
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





def _get_conserved_aa_exact_symbol(aa_dict, j_symbol):
    if "J-PHE" in aa_dict[j_symbol] and aa_dict[j_symbol]["J-PHE"] == "F":
        return "F"
    elif "J-TRP" in aa_dict[j_symbol] and aa_dict[j_symbol]["J-TRP"] == "W":
        return "W"
    elif "J-CYS" in aa_dict[j_symbol] and aa_dict[j_symbol]["J-CYS"] == "C":
        return "C"
    elif "J-VAL" in aa_dict[j_symbol] and aa_dict[j_symbol]["J-VAL"] == "V":
        return "V"



def is_valid_extension(original, extension):
    if extension.startswith(original): # todo what about AV/DV genes with no -1?
        if original[-1].isnumeric() and extension[len(original)].isnumeric():
            return False # cannot concatenate numbers if original ends with a number (e.g., TRAV1 to TRAV13 is not valid)
        if original[-1].isnumeric() and extension[len(original)] == "P":
            return False # do not add 'pseudogene'
        return True
    return False

def _get_all_extended_symbols(j_symbol, aa_dict, gene=""):
    return sorted(list({key for key in aa_dict.keys() if is_valid_extension(j_symbol, key) if gene in key}))


def _resolve_conserved_aa_from_partial_j_symbol(j_symbol, aa_dict, log_failures):
    extended_symbols = _get_all_extended_symbols(j_symbol, aa_dict, gene="J")
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


def _lookup_conserved_aa(j_symbol, aa_dict, log_failures):
    if j_symbol in aa_dict:
        return _get_conserved_aa_exact_symbol(aa_dict, j_symbol)

    return _resolve_conserved_aa_from_partial_j_symbol(j_symbol, aa_dict, log_failures)


def _get_aa_dict(gene_type, species):
    assert gene_type[0:2] in ("TR", "IG")
    species = species.lower().replace(" ", "")

    if gene_type[0:2] == "TR" and species == "homosapiens":
        return HOMOSAPIENS_TR_AA_SEQUENCES
    elif gene_type[0:2] == "IG" and species == "homosapiens":
        return HOMOSAPIENS_IG_AA_SEQUENCES
    elif gene_type[0:2] == "TR" and species == "musmusculus":
        return MUSMUSCULUS_TR_AA_SEQUENCES


def lon_get_conserved_aa_for_j_symbol_for_species(j_symbol, species, log_failures):
    if j_symbol.startswith("TR") or j_symbol.startswith("IG"):
        aa_dict = _get_aa_dict(j_symbol, species)

        if aa_dict is not None:
            return _lookup_conserved_aa(j_symbol, aa_dict, log_failures)

        else:
            raise ValueError(f"Failed to determine the conserved trailing amino acid for J symbol {j_symbol}: this feature "
                f"is not supported for {j_symbol[0:2]} genes of species {species}.")
    else:
        raise ValueError(f"Failed to determine the conserved trailing amino acid for J symbol {j_symbol}: symbol is not formatted correctly."
                         f"Please use tt.tr.standardize or tt.ig.standardize to correct the symbol.")

def get_conserved_aa_for_locus_for_species(locus, species, log_failures):
    if locus is not None:
        return lon_get_conserved_aa_for_j_symbol_for_species(locus, species, log_failures)

def get_conserved_aa(j_symbol, locus, species, log_failures):
    conserved_aa = lon_get_conserved_aa_for_j_symbol_for_species(j_symbol=j_symbol, species=species, log_failures=log_failures)  # returns None if aa is ambiguous

    if conserved_aa is None:
        conserved_aa = get_conserved_aa_for_locus_for_species(locus=locus, species=species, log_failures=log_failures)

    return conserved_aa


def get_all_aa_seqs_for_symbol(symbol, species, gene=""):
    aa_dict = _get_aa_dict(symbol, species)

    if aa_dict is not None:
        if symbol in aa_dict:
            return {symbol: aa_dict[symbol]}

        extended_symbols = _get_all_extended_symbols(symbol, aa_dict, gene)
        return {ext_symbol: aa_dict[ext_symbol] for ext_symbol in extended_symbols}

    return dict()


def _get_genes_to_alleles(alleles):
    genes_to_alleles = dict()

    for allele in alleles:
        gene = remove_allele(allele)

        if gene not in genes_to_alleles:
            genes_to_alleles[gene] = [allele]
        else:
            genes_to_alleles[gene].append(allele)

    return genes_to_alleles


def collapse_aas_per_gene(symbol_to_aa_dict):
    genes_to_alleles = _get_genes_to_alleles(symbol_to_aa_dict.keys())
    genes_to_aa = dict()

    for gene in genes_to_alleles:
        all_genes_same_aa = True
        gene_aa_dict = None

        for allele in genes_to_alleles[gene]:
            if gene_aa_dict is None:
                gene_aa_dict = symbol_to_aa_dict[allele]
            else:
                if gene_aa_dict != symbol_to_aa_dict[allele]:
                    all_genes_same_aa = False

        if all_genes_same_aa:
            genes_to_aa[gene] = gene_aa_dict
        else:
            for allele in genes_to_alleles[gene]:
                genes_to_aa[allele] = symbol_to_aa_dict[allele]

    return genes_to_aa


def get_all_collapsed_aa_seqs_for_symbol(symbol, species, gene=""):
    all_aas = get_all_aa_seqs_for_symbol(symbol, species, gene=gene)

    if "*" in symbol:
        return all_aas

    return collapse_aas_per_gene(all_aas)
