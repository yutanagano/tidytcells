import re
from tidytcells._utils.conserved_aa_lookup import is_valid_extension, get_all_aa_seqs_for_symbol


def score_j_alignment(seq, j_region, j_offset):
    penalty = -1.5
    max_mismatches = 1

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

    if any(matching_aas):
        # any mismatches at the beginning of the J are considered random deletion, no penalty
        matching_aas = matching_aas[matching_aas.index(True):]

        while matching_aas.count(False) > max_mismatches:
            matching_aas = matching_aas[matching_aas.index(False)+1:]

            if len(matching_aas) > 0 and any(matching_aas):
                matching_aas = matching_aas[matching_aas.index(True):]
            else:
                return -1

        # if matching_aas.count(False) > 1:
        #     return -1

        score = sum([1 if match == True else penalty for match in matching_aas])

        if score > 1 and int(score) > score:
            pass
        if score > 3:
            pass

        return -1 if score < 0 else score

    return -1

def valid_j_anchor(seq, j_region, j_conserved_idx, j_offset):
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


def get_best_j_alignment(seq, j_region, j_conserved_idx=None):
    best_offset = None
    best_score = -1

    for j_offset in range(len(seq) - len(j_region), len(seq)):
        if j_conserved_idx is None or valid_j_anchor(seq, j_region, j_conserved_idx, j_offset):
            score = score_j_alignment(seq, j_region, j_offset)

            if score > best_score:
                best_score = score
                best_offset = j_offset

    return best_offset


def print_j_alignment(seq, j_region, j_offset):
    if j_offset is None:
        print("No alignment:")
        print(seq)
        print(j_region)
        print("-------------------------------")
    elif j_offset >= 0:
        print(seq)
        print((" " * j_offset) + j_region)
        print("-------------------------------")
    else:
        print((" " * abs(j_offset)) + seq)
        print(j_region)
        print("-------------------------------")



def align_to_all_j(seq, j_symbol, species, min_score=3):
    aa_dict = get_all_aa_seqs_for_symbol(j_symbol, species, gene="J")

    best_genes = []
    best_score = min_score

    for gene, gene_seqs_dict in aa_dict.items():
        j_region = gene_seqs_dict["J-REGION"]

        if "J-MOTIF" in gene_seqs_dict:
            j_offset = get_best_j_alignment(seq, j_region, j_region.index(gene_seqs_dict["J-MOTIF"]))

            if j_offset is not None:
                score = score_j_alignment(seq, j_region, j_offset)

                if score > best_score:
                    best_score = score
                    best_genes = [gene]
                elif score == best_score:
                    best_genes.append(gene)

    best_score = best_score if len(best_genes) > 0 else -1

    return best_genes, best_score


def align_j(seq, j_symbol, species):
    aa_dict = get_all_aa_seqs_for_symbol(j_symbol, species)

    for gene, gene_seqs_dict in aa_dict.items():
        j_region = gene_seqs_dict["J-REGION"]
        j_offset = get_best_j_alignment(seq, j_region, j_region.index(gene_seqs_dict["J-MOTIF"]))

        if j_offset is None:
            print(f"Could not align with: {gene}       {seq}   {j_region}")
            align_to_all_j(seq, j_symbol, species)
            print("-----")
            pass





def score_v_alignment(seq, v_region, v_offset):
    assert v_offset < 0, "V offset can only be smaller than 0"
    penalty = -1.5

    v_region_align = v_region[v_offset:]
    seq_align = seq[:len(v_region_align)]

    matching_aas = [a == b for a, b in zip(seq_align, v_region_align)]

    if any(matching_aas):
        matching_aas = matching_aas[:len(matching_aas) - matching_aas[::-1].index(True)]
        return sum([1 if m else penalty for m in matching_aas])

    return -1


def valid_v_anchor(seq, v_region, v_conserved_idx_neg, v_offset):
    assert v_conserved_idx_neg < 0
    assert v_offset < 0

    # case 1: cdr3 seq starts right after conserved C
    if v_offset - 1 == v_conserved_idx_neg:
        return True

    # case 2: conserved C is somewhere in the seq, make sure seq head matches V
    if v_offset <= v_conserved_idx_neg:
        seq_conserved = seq[0:(v_conserved_idx_neg - v_offset) + 1]
        v_conserved = v_region[v_conserved_idx_neg - len(seq_conserved) + 1:v_conserved_idx_neg + 1]

        return seq_conserved == v_conserved

    return False


def get_best_v_alignment(seq, v_region, v_conserved_idx_neg):
    best_offset = None
    best_score = -1

    for v_offset in range(-1, v_conserved_idx_neg - len(seq), -1):
        if v_conserved_idx_neg is None or valid_v_anchor(seq, v_region, v_conserved_idx_neg, v_offset):
            score = score_v_alignment(seq, v_region, v_offset)
            if score > best_score:
                best_score = score
                best_offset = v_offset

    return best_offset


def align_to_all_v(seq, v_symbol, species, min_score=2):
    aa_dict = get_all_aa_seqs_for_symbol(v_symbol, species, gene="V")

    best_genes = []
    best_score = min_score

    for gene, gene_seqs_dict in aa_dict.items():
        v_region = gene_seqs_dict["V-REGION"]
        fr3 = gene_seqs_dict["FR3-IMGT"]

        v_conserved_idx_neg = v_region.index(fr3) + len(fr3) - (len(v_region) + 1)

        v_offset = get_best_v_alignment(seq, v_region, v_conserved_idx_neg)

        if v_offset is not None:
            score = score_v_alignment(seq, v_region, v_offset)

            if score > best_score:
                best_score = score
                best_genes = [gene]
            elif score == best_score:
                best_genes.append(gene)

            # if score >= best_score:
            #     print("score", best_score)
            #     print_v_alignment(seq, v_region, v_offset)

    best_score = best_score if len(best_genes) > 0 else -1

    return best_genes, best_score


def print_v_alignment(seq, v_region, v_offset):
    print(seq)
    print(v_region[v_offset:])
    print("--------------")
#
#
# def score_v_alignment(seq, v_region, v_offset):
#     penalty = -1
#
#     if v_offset >= 0:
#         v_start = v_offset
#     else:
#         v_start = len(v_region) + v_offset
#
#     alignment_length = min(len(seq), max(0, len(v_region) - v_start))
#
#     if alignment_length <= 0:
#         return -1
#
#     seq_align = seq[:alignment_length]
#     v_region_align = v_region[v_start:v_start + alignment_length]
#
#     matching_aas = [seq_aa == v_aa for seq_aa, v_aa in zip(seq_align, v_region_align)]
#
#     if any(matching_aas):
#         matching_aas = matching_aas[matching_aas.index(True):]
#         score = sum([1 if match == True else penalty for match in matching_aas])
#         return score
#
#     return -1
#
#
#
#
#
# def get_best_v_alignment(seq, v_region, v_conserved_idx=None):
#     best_offset = None
#     best_score = -1
#
#     for v_offset in range(-len(v_region), 1):
#         if v_conserved_idx is None or valid_v_anchor(seq, v_region, v_conserved_idx, v_offset):
#             score = score_v_alignment(seq, v_region, v_offset)
#
#             if score > best_score:
#                 best_score = score
#                 best_offset = v_offset
#
#     return best_offset
#
#
# def print_v_alignment(seq, v_region, v_offset):
#     if v_offset is None:
#         print("No alignment:")
#         print(seq)
#         print(v_region)
#         print("-------------------------------")
#     else:
#         if v_offset >= 0:
#             v_start = v_offset
#         else:
#             v_start = len(v_region) + v_offset
#
#         if v_start < 0:
#             print((" " * abs(v_start)) + v_region)
#             print(seq)
#             print("-------------------------------")
#         else:
#             print(v_region)
#             print((" " * v_start) + seq)
#             print("-------------------------------")
