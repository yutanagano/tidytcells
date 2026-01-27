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



def get_is_valid_locus_gene_fn(locus: str, gene: str):
    '''
    Utility funciton to check if a given gene is valid for a given locus, considering the shared TRAV/DV genes, and
    the fact that TRA chains frequently use TRDV genes (even if they're not called TRAV/DV)
    source: doi.org/10.1016/j.immuno.2025.100058
    '''
    if gene == "V":
        if locus == "TRA":
            return lambda x: x.startswith("TRAV") or x.startswith("TRDV")

        if locus == "TRD":
            return lambda x: x.startswith("TRDV") or (x.startswith("TRAV") and "/DV" in x)

    if len(locus) == 3:
        return lambda x: x.startswith(locus + gene)

    return lambda x: x.startswith(locus) and gene in x

def get_compatible_symbols(symbol, aa_dict, gene, locus, enforce_functional):
    '''
    For a given gene symbol, return all alleles which are considered to be valid extensions of the gene symbol.
    A gene symbol can be an allele (return itself), gene, subgroup, or even locus (TRA) or receptor (TR)
    '''
    is_valid_locus_gene = get_is_valid_locus_gene_fn(locus, gene)

    return [
        candidate
        for candidate in aa_dict.keys()
        if (not enforce_functional or aa_dict[candidate]["functionality"] in {"F", "[F]", "(F)"})
        and is_valid_locus_gene(candidate)
        and candidate.startswith(symbol)
        and not (symbol[-1].isnumeric() and candidate[len(symbol)].isnumeric())
        and not (symbol[-1].isnumeric() and candidate[len(symbol)] == "P")  # do not add 'pseudogene'
    ]

def _get_genes_to_alleles(alleles):
    '''Create mapping from gene to all its alleles'''
    genes_to_alleles = dict()

    for allele in alleles:
        gene = allele.rsplit("*")[0]

        if gene not in genes_to_alleles:
            genes_to_alleles[gene] = [allele]
        else:
            genes_to_alleles[gene].append(allele)

    return genes_to_alleles


def collapse_aa_dict_per_gene(symbol_to_aa_dict):
    '''
    Reduce the number of items in the sequence_aa dict by collapsing all alleles at the gene level
    for any gene where all alleles have the same sequence (keep separate alleles if sequences differ)
    '''
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


def _get_max_score_from_matches(matching_aas, max_mismatches, mismatch_penalty, gene):
    '''
    Utility function for computing the alignment score based on a boolean list of matching/non-matching amino acids.

    Alignment score = number of matches
    - Ignore any mismatches at the beginning of J or end of V (considered lost in V(D)J recombination);
      iteratively trim the boolean list until the remaining number of mismatches is under the max_mismatches cutoff
    - Subtract a penalty for remaining mismatches
    '''
    matching_aas = matching_aas[::-1] if gene == "V" else matching_aas

    if any(matching_aas):
        # any mismatches at the beginning of the J or end of the V are considered random deletion, no penalty
        matching_aas = matching_aas[matching_aas.index(True):]

        while matching_aas.count(False) > max_mismatches:
            matching_aas = matching_aas[matching_aas.index(False) + 1:]

            if len(matching_aas) > 0 and any(matching_aas):
                matching_aas = matching_aas[matching_aas.index(True):]
            else:
                return -1

        score = sum([1 if match == True else mismatch_penalty for match in matching_aas])

        if score > 0:
            return score

    return -1


def get_j_alignment_score(seq, j_region, j_offset, max_mismatches = 1, mismatch_penalty = -1.5):
    '''
    Determines the alignment score given a J region, CDR3 sequence and alignment offset.
    - Slice out the overlapping part of the alignment between the CDR3 and J region
    - For each amino acid in the slice, determine if they match (True/False)
    - Use _get_max_score_from_matches to compute the alignment score from the list of boolean matches
    '''

    # figure out start positions
    s_start = max(j_offset, 0)
    j_start = max(-j_offset, 0)

    # truncate to same length
    alignment_length = min(len(seq) - s_start, len(j_region) - j_start)

    if alignment_length <= 0:
        return -1

    seq_align = seq[s_start:s_start + alignment_length]
    j_region_align = j_region[j_start:j_start + alignment_length]

    matching_aas = [seq_aa == j_aa for seq_aa, j_aa in zip(seq_align, j_region_align)]

    return _get_max_score_from_matches(matching_aas, max_mismatches, mismatch_penalty, gene="J")


def get_v_alignment_score(seq, v_region, v_offset, max_mismatches = 1, mismatch_penalty = -1.5):
    '''
    Determines the alignment score given a V region, CDR3 sequence and alignment offset.
    - Slice out the overlapping part of the alignment between the CDR3 and V region
    - For each amino acid in the slice, determine if they match (True/False)
    - Use _get_max_score_from_matches to compute the alignment score from the list of boolean matches
    '''

    assert v_offset < 0, "V offset can only be smaller than 0"

    v_region_align = v_region[v_offset:]
    seq_align = seq[:len(v_region_align)]

    matching_aas = [a == b for a, b in zip(seq_align, v_region_align)]

    return _get_max_score_from_matches(matching_aas, max_mismatches, mismatch_penalty, gene="V")



def j_anchor_is_valid(seq, j_region, j_conserved_idx, j_offset):
    '''
    Given a J alignment offset, this function checks if the CDR3 sequence does not contradict the conserved anchor
    residues (anything after the conserved F/W).
    '''

    # safety check: index exists in J region
    if not (0 <= j_conserved_idx < len(j_region)):
        return False

    seq_anchor_start = j_offset + j_conserved_idx

    # case 1: conserved anchor would be just after the last aa of sequence
    if seq_anchor_start == len(seq):
        return True

    # case 2: conserved anchor is within seq, check if tail matches J region
    seq_tail = seq[seq_anchor_start:]
    j_tail = j_region[j_conserved_idx:j_conserved_idx + len(seq_tail)]
    return seq_tail == j_tail


def v_anchor_is_valid(seq, v_region, v_conserved_idx_neg, v_offset):
    '''
    Given a V alignment offset, this function checks if the CDR3 sequence does not contradict the conserved anchor
    residues (anything before the conserved C).
    '''
    assert v_conserved_idx_neg < 0
    assert v_offset < 0

    # case 1: conserved C is outside CDR3 seq (reconstruction)
    if v_offset - 1 >= v_conserved_idx_neg:
        return True

    # case 2: conserved C is somewhere in the seq, make sure seq head matches V
    seq_conserved = seq[0:(v_conserved_idx_neg - v_offset) + 1]
    v_conserved = v_region[v_conserved_idx_neg - len(seq_conserved) + 1:v_conserved_idx_neg + 1]

    return seq_conserved == v_conserved


def get_best_j_alignment_for_region(seq, j_region, j_conserved_idx, max_mismatches = 1, mismatch_penalty = -1.5):
    '''
    Get the best alignment index for one given J region

    :param seq: CDR3 sequence to be aligned
    :param v_region: J region sequence
    :param v_conserved_idx_neg: index of the conserved F/W in the J region
    :param max_mismatches: maximum number of mismatches allowed
    :param mismatch_penalty: penalty subtracted for single mismatch
    :return: alignment offset of the best-scoring alignment, plus the score
    '''
    best_offset = None
    best_score = -1

    for j_offset in range(len(seq) - len(j_region), len(seq)):
        if j_anchor_is_valid(seq, j_region, j_conserved_idx, j_offset):
            score = get_j_alignment_score(seq, j_region, j_offset, max_mismatches=max_mismatches, mismatch_penalty=mismatch_penalty)

            if score > best_score:
                best_score = score
                best_offset = j_offset

    return best_offset, best_score



def get_best_v_alignment_for_region(seq, v_region, v_conserved_idx_neg, max_mismatches = 1, mismatch_penalty = -1.5):
    '''
    Get the best alignment index for one given V region

    :param seq: CDR3 sequence to be aligned
    :param v_region: V region sequence
    :param v_conserved_idx_neg: index of the conserved C in the V region expressed as negative number
                                (number of aa's from the right)
    :param max_mismatches: maximum number of mismatches allowed
    :param mismatch_penalty: penalty subtracted for single mismatch
    :return: alignment offset of the best-scoring alignment, plus the score
    '''
    best_offset = None
    best_score = -1

    for v_offset in range(-1, v_conserved_idx_neg - len(seq), -1):
        if v_anchor_is_valid(seq, v_region, v_conserved_idx_neg, v_offset):
            score = get_v_alignment_score(seq, v_region, v_offset, max_mismatches=max_mismatches, mismatch_penalty=mismatch_penalty)

            if score > best_score:
                best_score = score
                best_offset = v_offset

    return best_offset, best_score


def align_j_regions(seq, j_aa_dict, min_j_score, mismatch_penalty, max_mismatches):
    '''
    Get all best alignments across all J genes

    :param seq: CDR3 sequence to be aligned
    :param j_aa_dict: dictionary of J gene sequence details (J-MOTIF + J-REGION) per allele/gene
    :param min_j_score: minimum overlap score required
    :param mismatch_penalty: penalty subtracted for single mismatch
    :param max_mismatches: maximum number of mismatches allowed
    :return: alignment details for all J genes with the best alignment scores
             this includes 'offset' indices (difference between conserved aa idx and end of aligned sequence idx)
             needed to calculate the corrected sequence
    '''
    best_score = min_j_score
    best_alignments = []

    for gene, gene_seqs_dict in j_aa_dict.items():
        j_region = gene_seqs_dict["J-REGION"]

        if "J-MOTIF" in gene_seqs_dict:
            j_conserved_idx = j_region.index(gene_seqs_dict["J-MOTIF"])

            j_offset, score = get_best_j_alignment_for_region(seq, j_region, j_conserved_idx, max_mismatches, mismatch_penalty)
            cur_alinment = {"gene": gene, "j_region": j_region, "j_conserved_idx": j_conserved_idx, "j_offset": j_offset,
                            "score": score}

            if score > best_score:
                best_score = score
                best_alignments = [cur_alinment]
            elif score == best_score:
                best_alignments.append(cur_alinment)

    return best_alignments


def get_conserved_c_idx(gene, v_gene_seqs_dict):
    '''
    Compute the conserved 'C' index for a given V gene, expressed as negative number (number of aa's from the right)
    '''
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


def align_v_regions(seq, v_aa_dict, min_v_score, mismatch_penalty, max_mismatches):
    '''
    Get all best alignments across all V genes

    :param seq: CDR3 sequence to be aligned
    :param v_aa_dict: dictionary of V gene sequence details (V-REGION) per allele/gene
    :param min_v_score: minimum overlap score required
    :param mismatch_penalty: penalty subtracted for single mismatch
    :param max_mismatches: maximum number of mismatches allowed
    :return: alignment details for all V genes with the best alignment scores
             this includes 'offset' indices (difference between V region end and CDR3 start)
             needed to calculate the corrected sequence
    '''
    best_score = min_v_score
    best_alignments = []

    for gene, gene_seqs_dict in v_aa_dict.items():
        v_conserved_idx_neg = get_conserved_c_idx(gene, gene_seqs_dict)

        if v_conserved_idx_neg is not None:
            v_offset, score = get_best_v_alignment_for_region(seq, gene_seqs_dict["V-REGION"], v_conserved_idx_neg, max_mismatches, mismatch_penalty)
            cur_alinment = {"gene": gene, "v_region": gene_seqs_dict["V-REGION"], "v_conserved_idx": v_conserved_idx_neg, "v_offset": v_offset,
                            "score": score}

            if score > best_score:
                best_score = score
                best_alignments = [cur_alinment]
            elif score == best_score:
                best_alignments.append(cur_alinment)

    return best_alignments

