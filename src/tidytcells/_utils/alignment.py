import logging
from collections.abc import Callable, Iterable

logger = logging.getLogger(__name__)


def get_is_valid_locus_gene_fn(locus: str, gene: str) -> Callable[[str], bool]:
    """
    Utility funciton to check if a given gene is valid for a given locus, considering the shared TRAV/DV genes, and
    the fact that TRA chains frequently use TRDV genes (even if they're not called TRAV/DV)
    source: doi.org/10.1016/j.immuno.2025.100058
    """
    if gene == "V":
        if locus == "TRA" or locus == "TRA/D":
            return lambda symbol: symbol.startswith("TRAV") or symbol.startswith("TRDV")

        if locus == "TRD":
            return lambda symbol: (
                symbol.startswith("TRDV")
                or (symbol.startswith("TRAV") and "/DV" in symbol)
            )

    if len(locus) == 3:
        return lambda symbol: symbol.startswith(locus + gene)

    return lambda symbol: symbol.startswith(locus) and gene in symbol


def get_is_valid_locus_multi_gene_fn(locus: str, gene: str) -> Callable[[str], bool]:
    """
    wrapper around get_is_valid_locus_gene_fn for multi-gene strings
    """

    is_valid_locus_gene_fn = get_is_valid_locus_gene_fn(locus, gene)
    return lambda symbol: all(
        is_valid_locus_gene_fn(sub_symbol.strip()) for sub_symbol in symbol.split(",")
    )


def is_valid_functionality(
    aa_dict: dict[str, dict[str, str]], candidate: str, enforce_functional: bool
) -> bool:
    if (not enforce_functional) or (
        aa_dict[candidate]["functionality"] in {"F", "[F]", "(F)"}
    ):
        return True

    return False


def get_compatible_symbols(
    symbol: str,
    aa_dict: dict[str, dict[str, str]],
    gene: str,
    locus: str,
    enforce_functional: bool,
) -> list[str]:
    """
    For a given gene symbol, return all alleles which are considered to be valid extensions of the gene symbol.
    A gene symbol can be an allele (return itself), gene, subgroup, or even locus (TRA) or receptor (TR)
    """
    if symbol in aa_dict.keys():
        if is_valid_functionality(aa_dict, symbol, enforce_functional):
            return [symbol]
        return []

    is_valid_locus_gene = get_is_valid_locus_gene_fn(locus, gene)

    return [
        candidate
        for candidate in aa_dict.keys()
        if is_valid_functionality(aa_dict, candidate, enforce_functional)
        and is_valid_locus_gene(candidate)
        and candidate.startswith(symbol)
        and not (symbol[-1].isnumeric() and candidate[len(symbol)].isnumeric())
        and not (
            symbol[-1].isnumeric() and candidate[len(symbol)] == "P"
        )  # do not add 'pseudogene'
    ]


def _get_genes_to_alleles(alleles: Iterable[str]) -> dict[str, list[str]]:
    """Create mapping from gene to all its alleles"""
    genes_to_alleles = dict()

    for allele in alleles:
        gene = allele.rsplit("*")[0]

        if gene not in genes_to_alleles:
            genes_to_alleles[gene] = [allele]
        else:
            genes_to_alleles[gene].append(allele)

    return genes_to_alleles


def collapse_aa_dict_per_gene(
    symbol_to_aa_dict: dict[str, dict[str, str]],
) -> dict[str, dict[str, str]]:
    """
    Reduce the number of items in the sequence_aa dict by collapsing all alleles at the gene level
    for any gene where all alleles have the same sequence (keep separate alleles if sequences differ)
    """
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


def _get_max_score_from_matches(
    matching_aas: list[bool],
    max_mismatches: int,
    mismatch_penalty: float,
    gene: str,
) -> float:
    """
    Utility function for computing the alignment score based on a boolean list of matching/non-matching amino acids.

    Alignment score = number of matches
    - Ignore any mismatches at the beginning of J or end of V (considered lost in V(D)J recombination);
      iteratively trim the boolean list until the remaining number of mismatches is under the max_mismatches cutoff
    - Subtract a penalty for remaining mismatches
    """
    matching_aas = matching_aas[::-1] if gene == "V" else matching_aas

    if any(matching_aas):
        # any mismatches at the beginning of the J or end of the V are considered random deletion, no penalty
        matching_aas = matching_aas[matching_aas.index(True) :]

        while matching_aas.count(False) > max_mismatches:
            matching_aas = matching_aas[matching_aas.index(False) + 1 :]

            if len(matching_aas) > 0 and any(matching_aas):
                matching_aas = matching_aas[matching_aas.index(True) :]
            else:
                return -1

        score = sum([1 if match else mismatch_penalty for match in matching_aas])

        if score > 0:
            return score

    return -1


def get_j_alignment_score(
    seq: str,
    j_region: str,
    j_offset: int,
    max_mismatches: int = 1,
    mismatch_penalty: float = -1.5,
) -> float:
    """
    Determines the alignment score given a J region, CDR3 sequence and alignment offset.
    - Slice out the overlapping part of the alignment between the CDR3 and J region
    - For each amino acid in the slice, determine if they match (True/False)
    - Use _get_max_score_from_matches to compute the alignment score from the list of boolean matches
    """

    # figure out start positions
    s_start = max(j_offset, 0)
    j_start = max(-j_offset, 0)

    # truncate to same length
    alignment_length = min(len(seq) - s_start, len(j_region) - j_start)

    if alignment_length <= 0:
        return -1

    seq_align = seq[s_start : s_start + alignment_length]
    j_region_align = j_region[j_start : j_start + alignment_length]

    matching_aas = [seq_aa == j_aa for seq_aa, j_aa in zip(seq_align, j_region_align)]

    return _get_max_score_from_matches(
        matching_aas, max_mismatches, mismatch_penalty, gene="J"
    )


def get_v_alignment_score(
    seq: str,
    v_region: str,
    v_offset: int,
    max_mismatches: int = 1,
    mismatch_penalty: float = -1.5,
) -> float:
    """
    Determines the alignment score given a V region, CDR3 sequence and alignment offset.
    - Slice out the overlapping part of the alignment between the CDR3 and V region
    - For each amino acid in the slice, determine if they match (True/False)
    - Use _get_max_score_from_matches to compute the alignment score from the list of boolean matches
    """

    assert v_offset < 0, "V offset can only be smaller than 0"

    v_region_align = v_region[v_offset:]
    seq_align = seq[: len(v_region_align)]

    matching_aas = [a == b for a, b in zip(seq_align, v_region_align)]

    return _get_max_score_from_matches(
        matching_aas, max_mismatches, mismatch_penalty, gene="V"
    )


def j_anchor_is_valid(
    seq: str, j_region: str, j_conserved_idx: int, j_offset: int
) -> bool:
    """
    Given a J alignment offset, this function checks if the CDR3 sequence does not contradict the conserved anchor
    residues (anything after the conserved F/W).
    """

    # safety check: index exists in J region
    if not (0 <= j_conserved_idx < len(j_region)):
        return False

    seq_anchor_start = j_offset + j_conserved_idx

    # case 1: conserved anchor would be just after the last aa of sequence
    if seq_anchor_start == len(seq):
        return True

    # case 2: conserved anchor is within seq, check if tail matches J region
    seq_tail = seq[seq_anchor_start:]
    j_tail = j_region[j_conserved_idx : j_conserved_idx + len(seq_tail)]
    return seq_tail == j_tail


def v_anchor_is_valid(
    seq: str, v_region: str, v_conserved_idx_neg: int, v_offset: int
) -> bool:
    """
    Given a V alignment offset, this function checks if the CDR3 sequence does not contradict the conserved anchor
    residues (anything before the conserved C).
    """
    assert v_conserved_idx_neg < 0
    assert v_offset < 0

    # case 1: conserved C is outside CDR3 seq (reconstruction)
    if v_offset - 1 >= v_conserved_idx_neg:
        return True

    # case 2: conserved C is somewhere in the seq, make sure seq head matches V
    seq_conserved = seq[0 : (v_conserved_idx_neg - v_offset) + 1]
    v_conserved = v_region[
        v_conserved_idx_neg - len(seq_conserved) + 1 : v_conserved_idx_neg + 1
    ]

    return seq_conserved == v_conserved


def get_best_j_alignment_for_region(
    seq: str,
    j_region: str,
    j_conserved_idx: int,
    max_mismatches: int = 1,
    mismatch_penalty: float = -1.5,
) -> tuple[int | None, float]:
    """
    Get the best alignment index for one given J region.

    Parameters
    ----------
    seq
        CDR3 sequence to be aligned.
    j_region
        J region sequence.
    j_conserved_idx
        Index of the conserved F/W in the J region.
    max_mismatches
        Maximum number of mismatches allowed.
    mismatch_penalty
        Penalty subtracted for single mismatch.

    Returns
    -------
    tuple[int | None, float]
        Alignment offset of the best-scoring alignment, plus the score.
    """
    best_offset = None
    best_score = -1

    for j_offset in range(len(seq) - len(j_region), len(seq)):
        if j_anchor_is_valid(seq, j_region, j_conserved_idx, j_offset):
            score = get_j_alignment_score(
                seq,
                j_region,
                j_offset,
                max_mismatches=max_mismatches,
                mismatch_penalty=mismatch_penalty,
            )

            if score > best_score:
                best_score = score
                best_offset = j_offset

    return best_offset, best_score


def get_best_v_alignment_for_region(
    seq: str,
    v_region: str,
    v_conserved_idx_neg: int,
    max_mismatches: int = 1,
    mismatch_penalty: float = -1.5,
) -> tuple[int | None, float]:
    """
    Get the best alignment index for one given V region.

    Parameters
    ----------
    seq
        CDR3 sequence to be aligned.
    v_region
        V region sequence.
    v_conserved_idx_neg
        Index of the conserved C in the V region expressed as negative number
        (number of aa's from the right).
    max_mismatches
        Maximum number of mismatches allowed.
    mismatch_penalty
        Penalty subtracted for single mismatch.

    Returns
    -------
    tuple[int | None, float]
        Alignment offset of the best-scoring alignment, plus the score.
    """
    best_offset = None
    best_score = -1

    for v_offset in range(-1, v_conserved_idx_neg - len(seq), -1):
        if v_anchor_is_valid(seq, v_region, v_conserved_idx_neg, v_offset):
            score = get_v_alignment_score(
                seq,
                v_region,
                v_offset,
                max_mismatches=max_mismatches,
                mismatch_penalty=mismatch_penalty,
            )

            if score > best_score:
                best_score = score
                best_offset = v_offset

    return best_offset, best_score


def align_j_regions(
    seq: str,
    j_aa_dict: dict[str, dict[str, str]],
    min_j_score: float,
    mismatch_penalty: float,
    max_mismatches: int,
) -> list[dict[str, str | int | float]]:
    """
    Get all best alignments across all J genes.

    Parameters
    ----------
    seq
        CDR3 sequence to be aligned.
    j_aa_dict
        Dictionary of J gene sequence details (J-MOTIF + J-REGION) per allele/gene.
    min_j_score
        Minimum overlap score required.
    mismatch_penalty
        Penalty subtracted for single mismatch.
    max_mismatches
        Maximum number of mismatches allowed.

    Returns
    -------
    list[dict[str, str | int | float]]
        Alignment details for all J genes with the best alignment scores.
        This includes 'offset' indices (difference between conserved aa idx and end of aligned sequence idx)
        needed to calculate the corrected sequence.
    """
    best_score = min_j_score
    best_alignments = []

    for gene, gene_seqs_dict in j_aa_dict.items():
        j_region = gene_seqs_dict["J-REGION"]

        if "J-MOTIF" in gene_seqs_dict:
            j_conserved_idx = j_region.index(gene_seqs_dict["J-MOTIF"])

            j_offset, score = get_best_j_alignment_for_region(
                seq, j_region, j_conserved_idx, max_mismatches, mismatch_penalty
            )
            cur_alinment = {
                "gene": gene,
                "j_region": j_region,
                "j_conserved_idx": j_conserved_idx,
                "j_offset": j_offset,
                "score": score,
            }

            if score > best_score:
                best_score = score
                best_alignments = [cur_alinment]
            elif score == best_score:
                best_alignments.append(cur_alinment)

    return best_alignments


def get_conserved_c_idx(gene: str, v_gene_seqs_dict: dict[str, str]) -> int | None:
    """
    Compute the conserved 'C' index for a given V gene, expressed as negative number (number of aa's from the right)
    """
    v_region = v_gene_seqs_dict["V-REGION"]
    v_conserved_idx_neg = None

    if "FR3-IMGT" in v_gene_seqs_dict:
        fr3 = v_gene_seqs_dict["FR3-IMGT"]
        try:
            v_conserved_idx_neg = v_region.index(fr3) + len(fr3) - (len(v_region) + 1)
        except ValueError:
            pass

    if gene.startswith("IG") and len(v_region) > 20 and "YYC" in v_region[-20:]:
        v_conserved_idx_neg = 0 - (len(v_region.rsplit("YYC", 1)[1]) + 1)

    if v_conserved_idx_neg is not None and v_region[v_conserved_idx_neg] == "C":
        return v_conserved_idx_neg


def align_v_regions(
    seq: str,
    v_aa_dict: dict[str, dict[str, str]],
    min_v_score: float,
    mismatch_penalty: float,
    max_mismatches: int,
) -> list[dict[str, str | int | float]]:
    """
    Get all best alignments across all V genes.

    Parameters
    ----------
    seq
        CDR3 sequence to be aligned.
    v_aa_dict
        Dictionary of V gene sequence details (V-REGION) per allele/gene.
    min_v_score
        Minimum overlap score required.
    mismatch_penalty
        Penalty subtracted for single mismatch.
    max_mismatches
        Maximum number of mismatches allowed.

    Returns
    -------
    list[dict[str, str | int | float]]
        Alignment details for all V genes with the best alignment scores.
        This includes 'offset' indices (difference between V region end and CDR3 start)
        needed to calculate the corrected sequence.
    """
    best_score = min_v_score
    best_alignments = []

    for gene, gene_seqs_dict in v_aa_dict.items():
        v_conserved_idx_neg = get_conserved_c_idx(gene, gene_seqs_dict)

        if v_conserved_idx_neg is not None:
            v_offset, score = get_best_v_alignment_for_region(
                seq,
                gene_seqs_dict["V-REGION"],
                v_conserved_idx_neg,
                max_mismatches,
                mismatch_penalty,
            )
            cur_alinment = {
                "gene": gene,
                "v_region": gene_seqs_dict["V-REGION"],
                "v_conserved_idx": v_conserved_idx_neg,
                "v_offset": v_offset,
                "score": score,
            }

            if score > best_score:
                best_score = score
                best_alignments = [cur_alinment]
            elif score == best_score:
                best_alignments.append(cur_alinment)

    return best_alignments
