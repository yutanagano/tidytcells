import logging
from collections.abc import Callable, Iterable

logger = logging.getLogger(__name__)


def get_is_valid_locus_gene_fn(locus: str, gene: str) -> Callable[[str], bool]:
    """
    Check if a given gene is valid for a given locus, considering the shared TRAV/DV genes, and
    the fact that TRA chains frequently use TRDV genes (even if they're not called TRAV/DV).
    source: doi.org/10.1016/j.immuno.2025.100058

    Parameters
    ----------
    locus
        Locus to check against.
    gene
        Gene to check against.

    Returns
    -------
    Callable[[str], bool]
        A function that checks whether a given symbol is valid for the locus and gene.
    """
    if gene == "V":
        if locus == "TRA" or locus == "TRA/D":
            return lambda symbol: symbol.startswith(("TRAV", "TRDV"))

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
    Wrapper around `get_is_valid_locus_gene_fn` for multi-gene strings.

    Parameters
    ----------
    locus
        Locus to check against.
    gene
        Gene to check against.

    Returns
    -------
    Callable[[str], bool]
        A function that checks whether all comma-separated symbols in a string are valid for the locus and gene.
    """

    is_valid_locus_gene_fn = get_is_valid_locus_gene_fn(locus, gene)
    return lambda symbol: all(
        is_valid_locus_gene_fn(sub_symbol.strip()) for sub_symbol in symbol.split(",")
    )


def is_valid_functionality(allele_data: dict, enforce_functional: bool) -> bool:
    return bool(not enforce_functional or allele_data.get("functionality") == "F")


def has_sequence_data(allele_data: dict) -> bool:
    return any(field != "functionality" for field in allele_data)


def get_compatible_symbols(
    symbol: str,
    allele_tree: dict[str, dict[str, dict]],
    gene: str,
    locus: str,
    enforce_functional: bool,
) -> list[str]:
    """
    For a given gene symbol, return all alleles which are considered to be valid extensions of the gene symbol.
    A gene symbol can be an allele (return itself), gene, subgroup, or even locus (TRA) or receptor (TR).

    Parameters
    ----------
    symbol
        Gene symbol to expand.
    allele_tree
        Nested dictionary of genes, alleles and their sequence details.
        Only alleles with known sequence data are returned.
    gene
        Gene to check validity against.
    locus
        Locus to check validity against.
    enforce_functional
        If ``True``, disallows non-functional genes / alleles.

    Returns
    -------
    list[str]
        All alleles which are considered to be valid extensions of the gene symbol.
    """
    if "*" in symbol:
        gene_name, allele_designation = symbol.split("*")
        allele_data = allele_tree.get(gene_name, {}).get(allele_designation)
        if (
            allele_data is not None
            and has_sequence_data(allele_data)
            and is_valid_functionality(allele_data, enforce_functional)
        ):
            return [symbol]
        return []

    is_valid_locus_gene = get_is_valid_locus_gene_fn(locus, gene)
    compatible_symbols = []

    for gene_name, alleles in allele_tree.items():
        if not gene_name.startswith(symbol):
            continue

        if symbol and symbol[-1].isnumeric() and len(gene_name) > len(symbol):
            next_char = gene_name[len(symbol)]
            if next_char.isnumeric() or next_char == "P":
                continue  # do not match longer gene names or 'pseudogene'

        if not is_valid_locus_gene(gene_name):
            continue

        for allele_designation, allele_data in alleles.items():
            if has_sequence_data(allele_data) and is_valid_functionality(
                allele_data, enforce_functional
            ):
                compatible_symbols.append(f"{gene_name}*{allele_designation}")

    return compatible_symbols


def collapse_alleles_with_same_aa(
    symbol_to_aa_dict: dict[str, dict[str, str]],
) -> dict[str, dict[str, str]]:
    """
    Reduce the number of items in the sequence_aa dict by collapsing all alleles with identical
    sequence information, regardless of gene. Each group is keyed by its sorted allele names
    joined with ", " (a single allele keeps its own name).

    Parameters
    ----------
    symbol_to_aa_dict
        Dictionary of amino acid sequence details per allele.

    Returns
    -------
    dict[str, dict[str, str]]
        The reduced dictionary with one item per group of identical alleles.
    """
    alleles_per_aa = {}

    for allele, aa_dict in symbol_to_aa_dict.items():
        aa_key = tuple(sorted(aa_dict.items()))

        if aa_key not in alleles_per_aa:
            alleles_per_aa[aa_key] = [allele]
        else:
            alleles_per_aa[aa_key].append(allele)

    return {", ".join(sorted(alleles)): symbol_to_aa_dict[alleles[0]]
            for alleles in alleles_per_aa.values()}


def _get_max_score_from_matches(
    matching_aas: list[bool],
    allow_mismatches: bool,
    mismatch_penalty: float,
    gene: str,
) -> int:
    """
    Utility function for computing the alignment score based on a boolean list of matching/non-matching amino acids.

    Alignment score = number of matches:
    - Ignore any mismatches at the beginning of J or end of V (considered lost in V(D)J recombination);
      iteratively trim the boolean list until the remaining number of mismatches is under the max_mismatches cutoff
    - Subtract a penalty for remaining mismatches

    Parameters
    ----------
    matching_aas
        Boolean list of matching (True) / non-matching (False) amino acids.
    allow_mismatches
        Whether mismatches are allowed.
    mismatch_penalty
        Penalty subtracted for single mismatch.
    gene
        Gene being aligned (``'V'`` or ``'J'``), used to orient the list.

    Returns
    -------
    int
        The best alignment score (``-1`` if no valid alignment).
    """
    if not any(matching_aas):
        return -1

    # orient the list from the side lost in recombination towards the anchor side
    if gene == "V":
        matching_aas = matching_aas[::-1]

    best_score = -1
    score = 0

    # extend the alignment one residue at a time from the anchor side,
    # keeping the best score seen
    for match in reversed(matching_aas):
        if match:
            score += 1
        elif allow_mismatches:
            score += mismatch_penalty
        else:
            break

        best_score = max(best_score, score)

    return best_score if best_score >= 1 else -1


def get_j_alignment_score(
    seq: str,
    j_region: str,
    j_offset: int,
    allow_mismatches: bool = True,
    mismatch_penalty: float = -1.5,
) -> int:
    """
    Determine the alignment score given a J region, CDR3 sequence and alignment offset.

    - Slice out the overlapping part of the alignment between the CDR3 and J region
    - For each amino acid in the slice, determine if they match (True/False)
    - Use `_get_max_score_from_matches` to compute the alignment score from the list of boolean matches

    Parameters
    ----------
    seq
        CDR3 sequence to be aligned.
    j_region
        J region sequence.
    j_offset
        Alignment offset.
    allow_mismatches
        Whether mismatches are allowed.
        Defaults to ``True``.
    mismatch_penalty
        Penalty subtracted for single mismatch.
        Defaults to ``-1.5``.

    Returns
    -------
    int
        The alignment score.
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
        matching_aas, allow_mismatches, mismatch_penalty, gene="J"
    )


def get_v_alignment_score(
    seq: str,
    v_region: str,
    v_offset: int,
    allow_mismatches: bool = True,
    mismatch_penalty: float = -1.5,
) -> int:
    """
    Determine the alignment score given a V region, CDR3 sequence and alignment offset.

    - Slice out the overlapping part of the alignment between the CDR3 and V region
    - For each amino acid in the slice, determine if they match (True/False)
    - Use `_get_max_score_from_matches` to compute the alignment score from the list of boolean matches

    Parameters
    ----------
    seq
        CDR3 sequence to be aligned.
    v_region
        V region sequence.
    v_offset
        Alignment offset.
    allow_mismatches
        Whether mismatches are allowed.
        Defaults to ``True``.
    mismatch_penalty
        Penalty subtracted for single mismatch.
        Defaults to ``-1.5``.

    Returns
    -------
    int
        The alignment score.
    """

    assert v_offset < 0, "V offset can only be smaller than 0"

    v_region_align = v_region[v_offset:]
    seq_align = seq[: len(v_region_align)]

    matching_aas = [a == b for a, b in zip(seq_align, v_region_align)]

    return _get_max_score_from_matches(
        matching_aas, allow_mismatches, mismatch_penalty, gene="V"
    )


def j_anchor_is_valid(
    seq: str,
    j_region: str,
    j_conserved_idx: int,
    j_offset: int,
    allow_flank_mismatches: bool = True,
) -> bool:
    """
    Given a J alignment offset, check if the CDR3 sequence
    - Matches at the conserved anchor position (if this falls inside the alignment)
    - Sufficiently matches the J end region (match at least 2/3 of amino acids for long alignments)

    Parameters
    ----------
    seq
        CDR3 sequence to be aligned.
    j_region
        J region sequence.
    j_conserved_idx
        Index of the conserved F/W in the J region.
    j_offset
        Alignment offset.
    allow_flank_mismatches
        Whether mismatches are allowed in the flanking region.
        Defaults to ``True``.

    Returns
    -------
    bool
        Whether the anchor is valid.
    """

    # safety check: index exists in J region
    if not (0 <= j_conserved_idx < len(j_region)):
        return False

    seq_anchor_start = j_offset + j_conserved_idx

    # case 1: conserved anchor would be after the last aa of sequence
    if seq_anchor_start >= len(seq):
        return True

    # case 2: conserved anchor is within seq, check if tail sufficiently matches J region
    seq_tail = seq[seq_anchor_start:]
    j_tail = j_region[j_conserved_idx : j_conserved_idx + len(seq_tail)]
    matches_tail = [seq_aa == j_aa for seq_aa, j_aa in zip(seq_tail, j_tail)]

    # conserved anchor *must* match
    if not matches_tail[0]:
        return False

    if allow_flank_mismatches:
        # permit some amount of mismatches in J end to account for SHM/errors
        # at least 2/3 of amino acids must match
        return 3 * sum(matches_tail) >= 2 * len(matches_tail)
    else:
        return all(matches_tail)


def v_anchor_is_valid(
    seq: str,
    v_region: str,
    v_conserved_idx_neg: int,
    v_offset: int,
    allow_flank_mismatches: bool = True,
) -> bool:
    """
    Given a V alignment offset, check if the CDR3 sequence
    - Matches at the conserved anchor position (if this falls inside the alignment)
    - Sufficiently matches the V frw region (1 mismatch for short alignments, or
      at least 2/3 of amino acids for long alignments)

    Parameters
    ----------
    seq
        CDR3 sequence to be aligned.
    v_region
        V region sequence.
    v_conserved_idx_neg
        Index of the conserved C in the V region expressed as negative number
        (number of aa's from the right).
    v_offset
        Alignment offset.
    allow_flank_mismatches
        Whether mismatches are allowed in the flanking region.
        Defaults to ``True``.

    Returns
    -------
    bool
        Whether the anchor is valid.
    """
    assert v_conserved_idx_neg < 0
    assert v_offset < 0

    # case 1: conserved C is outside CDR3 seq (reconstruction)
    if v_offset - 1 >= v_conserved_idx_neg:
        return True

    # case 2: conserved C is somewhere in the seq, make sure seq head (sufficiently) matches V
    seq_frw = seq[0 : (v_conserved_idx_neg - v_offset) + 1]

    v_conserved_pos = len(v_region) + v_conserved_idx_neg
    v_frw = v_region[v_conserved_pos - len(seq_frw) + 1 : v_conserved_pos + 1]

    matches_frw = [seq_aa == v_aa for seq_aa, v_aa in zip(seq_frw, v_frw)]

    # conserved anchor *must* match
    if not matches_frw[-1]:
        return False

    # permit some amount of mismatches in frw to account for SHM/errors
    if allow_flank_mismatches:
        # if conserved anchor matches, and frw contains at most 1 additional aa (len 2),
        # then we say 'anchor is valid' regardless of if that second aa mismatches
        # note that this only applies to the V side: the J FGxG motif is far more conserved
        if len(matches_frw) <= 2:
            return True

        # for longer flanks at least 2/3 of amino acids must match
        return 3 * sum(matches_frw) >= 2 * len(matches_frw)
    else:
        return all(matches_frw)


def get_best_j_alignment_for_region(
    seq: str,
    j_region: str,
    j_conserved_idx: int,
    allow_mismatches: bool = True,
    mismatch_penalty: float = -1.5,
) -> tuple[list[int], int]:
    """
    Get all equally best alignment offsets for one given J region.

    Parameters
    ----------
    seq
        CDR3 sequence to be aligned.
    j_region
        J region sequence.
    j_conserved_idx
        Index of the conserved F/W in the J region.
    allow_mismatches
        Whether mismatches are allowed.
        Defaults to ``True``.
    mismatch_penalty
        Penalty subtracted for single mismatch.
        Defaults to ``-1.5``.

    Returns
    -------
    tuple[list[int], int]
        List of alignment offsets that share the best score (empty if no valid alignment), plus the score.
    """
    best_offsets = []
    best_score = -1

    min_offset_j_region = len(seq) - len(j_region)
    min_offset_anchor = -j_conserved_idx

    for j_offset in range(max(min_offset_j_region, min_offset_anchor), len(seq)):
        if j_anchor_is_valid(
            seq,
            j_region,
            j_conserved_idx,
            j_offset,
            allow_flank_mismatches=allow_mismatches,
        ):
            score = get_j_alignment_score(
                seq,
                j_region,
                j_offset,
                allow_mismatches=allow_mismatches,
                mismatch_penalty=mismatch_penalty,
            )

            if score > best_score:
                best_score = score
                best_offsets = [j_offset]
            elif score == best_score and score > 0:
                best_offsets.append(j_offset)

    return best_offsets, best_score


def get_best_v_alignment_for_region(
    seq: str,
    v_region: str,
    v_conserved_idx_neg: int,
    allow_mismatches: bool = True,
    mismatch_penalty: float = -1.5,
) -> tuple[list[int], int]:
    """
    Get all equally best alignment offsets for one given V region.

    Parameters
    ----------
    seq
        CDR3 sequence to be aligned.
    v_region
        V region sequence.
    v_conserved_idx_neg
        Index of the conserved C in the V region expressed as negative number
        (number of aa's from the right).
    allow_mismatches
        Whether mismatches are allowed.
        Defaults to ``True``.
    mismatch_penalty
        Penalty subtracted for single mismatch.
        Defaults to ``-1.5``.

    Returns
    -------
    tuple[list[int], int]
        List of alignment offsets that share the best score (empty if no valid alignment), plus the score.
    """
    best_offsets = []
    best_score = -1

    for v_offset in range(v_conserved_idx_neg - len(seq) + 1, 0):
        if v_anchor_is_valid(
            seq,
            v_region,
            v_conserved_idx_neg,
            v_offset,
            allow_flank_mismatches=allow_mismatches,
        ):
            score = get_v_alignment_score(
                seq,
                v_region,
                v_offset,
                allow_mismatches=allow_mismatches,
                mismatch_penalty=mismatch_penalty,
            )

            if score > best_score:
                best_score = score
                best_offsets = [v_offset]
            elif score == best_score and score > 0:
                best_offsets.append(v_offset)

    return best_offsets, best_score


def align_j_regions(
    seq: str,
    j_aa_dict: dict[str, dict[str, str]],
    min_j_score: int,
    mismatch_penalty: float,
    allow_mismatches: bool,
) -> list[dict[str, str | int]]:
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
    allow_mismatches
        Whether mismatches are allowed.

    Returns
    -------
    list[dict[str, str | int]]
        Alignment details for all J genes with the best alignment scores, one entry per equally best offset.
        This includes 'offset' indices (difference between conserved aa idx and end of aligned sequence idx)
        needed to calculate the corrected sequence.
    """
    best_score = min_j_score
    best_alignments = []

    for gene, gene_seqs_dict in j_aa_dict.items():
        j_region = gene_seqs_dict["J-REGION"]

        if "J-MOTIF" in gene_seqs_dict:
            j_conserved_idx = j_region.index(gene_seqs_dict["J-MOTIF"])

            j_offsets, score = get_best_j_alignment_for_region(
                seq,
                j_region,
                j_conserved_idx,
                allow_mismatches=allow_mismatches,
                mismatch_penalty=mismatch_penalty,
            )
            cur_alignments = [
                {
                    "gene": gene,
                    "j_region": j_region,
                    "j_conserved_idx": j_conserved_idx,
                    "j_offset": j_offset,
                    "score": score,
                }
                for j_offset in j_offsets
            ]

            if score > best_score:
                best_score = score
                best_alignments = cur_alignments
            elif score == best_score:
                best_alignments.extend(cur_alignments)

    return best_alignments


def get_conserved_c_idx(gene: str, v_gene_seqs_dict: dict[str, str]) -> int | None:
    """
    Compute the conserved 'C' index for a given V gene, expressed as negative number (number of aa's from the right).

    Parameters
    ----------
    gene
        V gene symbol.
    v_gene_seqs_dict
        Dictionary of V gene sequence details.

    Returns
    -------
    int | None
        The conserved 'C' index, or ``None`` if it cannot be determined.
    """
    v_region = v_gene_seqs_dict["V-REGION-SHORT"]
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
    min_v_score: int,
    mismatch_penalty: float,
    allow_mismatches: bool,
) -> list[dict[str, str | int]]:
    """
    Get all best alignments across all V genes.

    Parameters
    ----------
    seq
        CDR3 sequence to be aligned.
    v_aa_dict
        Dictionary of V gene sequence details (V-REGION-SHORT) per allele/gene.
    min_v_score
        Minimum overlap score required.
    mismatch_penalty
        Penalty subtracted for single mismatch.
    allow_mismatches
        Whether mismatches are allowed.

    Returns
    -------
    list[dict[str, str | int]]
        Alignment details for all V genes with the best alignment scores.
        This includes 'offset' indices (difference between V region end and CDR3 start)
        needed to calculate the corrected sequence.
    """
    best_score = min_v_score
    best_alignments = []

    for gene, gene_seqs_dict in v_aa_dict.items():
        v_conserved_idx_neg = get_conserved_c_idx(gene, gene_seqs_dict)

        if v_conserved_idx_neg is not None:
            v_offsets, score = get_best_v_alignment_for_region(
                seq,
                gene_seqs_dict["V-REGION-SHORT"],
                v_conserved_idx_neg,
                allow_mismatches=allow_mismatches,
                mismatch_penalty=mismatch_penalty,
            )
            cur_alignments = [
                {
                    "gene": gene,
                    "v_region_short": gene_seqs_dict["V-REGION-SHORT"],
                    "v_conserved_idx": v_conserved_idx_neg,
                    "v_offset": v_offset,
                    "score": score,
                }
                for v_offset in v_offsets
            ]

            if score > best_score:
                best_score = score
                best_alignments = cur_alignments
            elif score == best_score:
                best_alignments.extend(cur_alignments)

    return best_alignments