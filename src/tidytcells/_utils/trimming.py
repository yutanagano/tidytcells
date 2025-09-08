import logging
import re

from tidytcells._resources import HOMOSAPIENS_TR_AA_SEQUENCES, HOMOSAPIENS_IG_AA_SEQUENCES, MUSMUSCULUS_TR_AA_SEQUENCES


logger = logging.getLogger(__name__)



def get_motifs_from_dict(aa_dict, motif_type, locus):
    locus = locus if locus is not None else ""
    return {value[motif_type] for key, value in aa_dict.items() if motif_type in value and key.startswith(locus)}

def get_motifs_for_tr_locus(tr_aa_dict, motif_type, locus):
    # TRAV/DV genes can be shared, recommended to add DV genes to AV (see https://doi.org/10.1016/j.immuno.2025.100058)
    # this is true for both mouse and human

    if locus in ("TRA", "TRD") and motif_type == "V-MOTIF":
        motifs = set()
        motifs.update(get_motifs_from_dict(tr_aa_dict, motif_type, locus="TRA"))
        motifs.update(get_motifs_from_dict(tr_aa_dict, motif_type, locus="TRD"))

        return motifs

    return get_motifs_from_dict(HOMOSAPIENS_TR_AA_SEQUENCES, motif_type, locus)

def get_conserved_motifs(locus, species, motif_type):
    assert motif_type in ("V-MOTIF", "J-MOTIF"), f"Unknown motif type: {motif_type}"

    if species == "homosapiens":
        if locus is None: # get all TR + IG motifs combined
            motifs = set()
            motifs.update(get_motifs_from_dict(HOMOSAPIENS_TR_AA_SEQUENCES, motif_type, locus))
            motifs.update(get_motifs_from_dict(HOMOSAPIENS_IG_AA_SEQUENCES, motif_type, locus))

            return motifs
        elif locus.startswith("TR"):
            return get_motifs_for_tr_locus(HOMOSAPIENS_TR_AA_SEQUENCES, motif_type, locus)
        elif locus.startswith("IG"):
            return get_motifs_from_dict(HOMOSAPIENS_IG_AA_SEQUENCES, motif_type, locus)
    elif species == "musmusculus":
        if locus is not None and locus.startswith("IG"):
            raise ValueError(f"Trimming is not supported for locus {locus} and species {species}.")

        return get_motifs_for_tr_locus(MUSMUSCULUS_TR_AA_SEQUENCES, motif_type, locus)
    else:
        raise ValueError(f"Trimming is not supported for species {species}.")


def trim_v_motif_start_recursive(orig_seq, v_motifs):
    v_motifs = {motif[1:] for motif in v_motifs}

    if len(v_motifs) == 0 or len(list(v_motifs)[0]) == 1:
        return orig_seq

    for motif in v_motifs:
        if orig_seq.startswith(motif):
            return orig_seq.replace(motif, "C", 1)

    return trim_v_motif_start_recursive(orig_seq, v_motifs)

def trim_j_motif_end_recursive(orig_seq, j_motifs):
    j_motifs = {motif[:-1] for motif in j_motifs}

    if len(j_motifs) == 0 or len(list(j_motifs)[0]) == 1:
        return orig_seq

    for motif in j_motifs:
        if orig_seq.endswith(motif):
            return orig_seq.rsplit(motif, 1)[0] + motif[-1]

    return trim_j_motif_end_recursive(orig_seq, j_motifs)


def trim_junction_start(orig_seq, locus, species):
    v_motifs = get_conserved_motifs(locus, species, "V-MOTIF")

    # case 1: few extra amino acids at start, specific known motifs
    #   trim these *only* if they appear at the start
    seq = trim_v_motif_start_recursive(orig_seq, v_motifs)

    if seq != orig_seq and seq.startswith("C"):
        return seq

    # case 2: long sequence provided (>4 extra aa's at start)
    #   look for conserved V motif of length 4 *anywhere* in the sequence. if found, trim
    #   due to the length, this rarely occurs randomly in CDR3
    for motif in v_motifs:
        if motif in orig_seq:
            new_seq = "C" + orig_seq.split(motif)[1]
            print("V trim orig_seq:", orig_seq, "new_seq:", new_seq)
            if len(new_seq) >= 6:
                print("--returned")
                return new_seq

    return orig_seq

def trim_junction_end(orig_seq, locus, species, conserved_aas):
    j_motifs = get_conserved_motifs(locus, species, "J-MOTIF")
    j_motifs = {motif for motif in j_motifs if motif[0] in conserved_aas}

    # case 1: few extra amino acids at start, specific known motifs
    #   trim these *only* if they appear at the start
    seq = trim_j_motif_end_recursive(orig_seq, j_motifs)

    if seq != orig_seq and seq[-1] in conserved_aas:
        return seq

    # case 2: long sequence provided (>4 extra aa's at start)
    #   look for conserved J motif of length 4 *anywhere* in the sequence. if found, trim
    #   due to the length, this rarely occurs randomly in CDR3
    # todo maybe have min length before this is applied?
    for motif in j_motifs:
        if motif in seq:
            new_seq = seq.split(motif)[0] + motif[0]
            print("J trim orig_seq:", orig_seq, "start_trimmed:", seq, "new_seq:", seq)
            if len(new_seq) >= 6:
                print("--returned")
                return new_seq

    return orig_seq

def get_expected_conserved_aas(junction_matching_regex):
    return junction_matching_regex.pattern.split("*")[1].strip("[]$")

def trim_junction(orig_seq, junction_matching_regex, locus, species):
    # do not trim if it is already a valid CDR3 junction, or if it contains no C (likely to only be CDR3 without conserved aa's)
    if junction_matching_regex.match(orig_seq) or "C" not in orig_seq:
        return orig_seq

    seq = trim_junction_start(orig_seq, locus, species)

    # if junction_matching_regex.match(seq): # WFGA
    #     return seq

    if seq.startswith("C"): # only continue end-trimming if start is 'complete'
        # todo can be made more efficient: first check from the right side if the position is in there
        conserved_aas = get_expected_conserved_aas(junction_matching_regex)
        seq = trim_junction_end(seq, locus, species, conserved_aas)

        if junction_matching_regex.match(seq):
            return seq

    # if any trimming failed, this is likely the CDR3 without conserved aa's; return original sequence
    return orig_seq

