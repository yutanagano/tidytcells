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
    # todo allow locus to be list? IGL/IGK?? or just use IG in that case?

    assert motif_type in ("V-MOTIF", "J-MOTIF", "V-CDR3-MOTIF", "J-CDR3-MOTIF"), f"Unknown motif type: {motif_type}"

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

def minimal_trim_junction_start(orig_seq, locus, species):
    # if junction / or CDR3
    if orig_seq[0] == "C" or startswith_cdr3_start_motif(orig_seq, locus, species):
        return orig_seq

    v_motifs = get_conserved_motifs(locus, species, "V-MOTIF")

    for i in range(1, len(orig_seq) - 5):
        if orig_seq[i] == "C":
            trimmed_seq = orig_seq[i:]

            motif_in_seq = orig_seq[:i+1]
            motif_in_seq = motif_in_seq[-4:] if len(motif_in_seq) >= 4 else motif_in_seq

            v_motifs_for_len = {motif[-len(motif_in_seq):] for motif in v_motifs}

            for known_motif in v_motifs_for_len:
                if known_motif == motif_in_seq:
                    return trimmed_seq

    return orig_seq

def minimal_trim_junction_end(orig_seq, locus, species, conserved_aas):
    # if junction / or CDR3
    if orig_seq[-1] in conserved_aas or endswith_cdr3_end_motif(orig_seq, locus, species):
        return orig_seq

    j_motifs = get_conserved_motifs(locus, species, "J-MOTIF")
    j_motifs = {motif for motif in j_motifs if motif[0] in conserved_aas}

    for i in range(2, len(orig_seq) - 4):
        if orig_seq[-i] in conserved_aas:
            trimmed_seq = orig_seq[:-i + 1]

            motif_in_seq = orig_seq[-i:]
            motif_in_seq = motif_in_seq[:4] if len(motif_in_seq) >= 4 else motif_in_seq

            j_motifs_for_len = {motif[:len(motif_in_seq)] for motif in j_motifs}

            for known_motif in j_motifs_for_len:
                if known_motif == motif_in_seq:
                    return trimmed_seq

    return orig_seq

def get_expected_conserved_aas(junction_matching_regex):
    return junction_matching_regex.pattern.split("[")[-1].rstrip("]$")


# def trim_junction(orig_seq, junction_matching_regex, conserved_aa, locus, species):
#     # do not trim if it is already a valid CDR3 junction, or if it contains no C (likely to only be CDR3 without conserved aa's)
#     if junction_matching_regex.match(orig_seq) or "C" not in orig_seq:
#         return orig_seq
#
#     seq = minimal_trim_junction_start(orig_seq, locus, species)
#
#     # if junction_matching_regex.match(seq): # WFGA
#     #     return seq
#
#     if seq.startswith("C"): # only continue end-trimming if start is 'complete'
#         conserved_aas = get_expected_conserved_aas(junction_matching_regex)
#
#         seq = minimal_trim_junction_end(seq, locus, species, conserved_aas)
#
#         if junction_matching_regex.match(seq):
#             # print(f"'{orig_seq}': '{seq}',")
#             return seq
#
#     # if any trimming failed, this is likely the CDR3 without conserved aa's; return original sequence
#     return orig_seq


def startswith_cdr3_start_motif(seq, locus, species):
    cdr3_start_motifs = get_conserved_motifs(locus, species, "V-CDR3-MOTIF")

    for motif in cdr3_start_motifs:
        if seq.startswith(motif):
            return True

    return False

def endswith_cdr3_end_motif(seq, locus, species):
    cdr3_end_motifs = get_conserved_motifs(locus, species, "J-CDR3-MOTIF")

    for motif in cdr3_end_motifs:
        if seq.endswith(motif):
            return True

    return False


def is_valid_junction(seq, junction_matching_regex, locus, species, check_motifs=True):
    if junction_matching_regex.match(seq):

        if check_motifs:
            return is_valid_cdr3(seq[1:-1], locus, species)
        else:
            return True

    return False


def is_valid_cdr3(seq, locus, species, min_len=4):
    if len(seq) >= min_len:
        if startswith_cdr3_start_motif(seq, locus, species) and endswith_cdr3_end_motif(seq, locus, species):
            return True

    return False


    # if junction_matching_regex.match(seq):
    #     if check_motifs:
    #         pass
    #     else:
    #         return True
    #
    #
    # return False


def process_junction(orig_seq, junction_matching_regex, conserved_aa, locus, species, trimming=True, check_motifs=True):
    # do not trim if it is already a valid CDR3 junction, or if it contains no C (likely to only be CDR3 without conserved aa's)
    # if junction_matching_regex.match(orig_seq): # or "C" not in orig_seq:
    #     # return orig_seq
    #     print("used to be return orig seq but this is not 'safe'?")
    #     pass

    # cases are checked in order of most likely occurrence + minimal interventions
    # case 1: is junction with leading/trailing, CAAAAF (already valid)
    # case 2: is CDR3 without leading/trailing, AAAA (reconstruct leading/trailing)
    # case 3: is junction with excess on 1 or both sides (fully recoverable substring) xxxCAAAAFxxx / CAAAAFxxx / xxxCAAAAF
    # case 4: is CDR3 with excess on 1 or both sides (mostly substring, reconstruct leading/trailing) AAAAFxxx / xxxCAAAA
    # case 5: is partial CDR3 / does not contain valid motifs (needs reconstruction?? only for TCR)

    # case 1: junction
    if is_valid_junction(orig_seq, junction_matching_regex, locus, species, check_motifs):
        print("valid junction:", orig_seq) # todo for the cdr3 motifs keep with conserved aa
        return orig_seq

    # case 2: cdr3
    if check_motifs and is_valid_cdr3(orig_seq, locus, species):
        # todo retrieve precise conserved aa if all_conserved_aas is ambiguous?, this one is just preferred
        print("valid cdr3 to junction:", orig_seq, ":", "C" + orig_seq + conserved_aa)

        return "C" + orig_seq + conserved_aa

    # trimming
    if trimming:
        seq = minimal_trim_junction_start(orig_seq, locus, species)
        all_conserved_aas = get_expected_conserved_aas(junction_matching_regex) # todo if all_conserved_aas >1, determine based on J motif
        seq = minimal_trim_junction_end(seq, locus, species, all_conserved_aas)
    else:
        seq = orig_seq

    # case 3
    if is_valid_junction(seq, junction_matching_regex, locus, species, check_motifs):
        print("valid junction after trimming:", orig_seq, ":", seq)
        return seq

    # case 4 (only add conserved aas to sides that were not trimmed)
    start_was_trimmed = not orig_seq.startswith(seq)
    end_was_trimmed = not orig_seq.endswith(seq)

    if startswith_cdr3_start_motif(seq, locus, species) and not start_was_trimmed:
        seq = "C" + seq

    if endswith_cdr3_end_motif(seq, locus, species) and not end_was_trimmed:
        seq + conserved_aa

    if is_valid_junction(seq, junction_matching_regex, locus, species, check_motifs):
        print("valid junction after trimming and restoring:", orig_seq, ":", seq)
        return seq

    # assert False, "incomplete or invalid CDR3?"
    # case 5: invalid CDR3

    print("incomplete or invalid CDR3:", orig_seq, f": {seq})" if seq != orig_seq else "")
    return None







    #
    # # 'CAELNAGNNRKLI' -> CAELNAGNNRKLIW if TCR
    #
    # if not start_was_trimmed and not end_was_trimmed:
    #     return orig_seq
    #
    # if start_was_trimmed or end_was_trimmed and junction_matching_regex.match(seq):
    #     return seq
    #
    # if not start_was_trimmed and end_was_trimmed and startswith_cdr3_start_motif(seq, locus, species):
    #     seq = "C" + seq
    #     if junction_matching_regex.match(seq):
    #         return seq
    #
    # if start_was_trimmed and not end_was_trimmed and endswith_cdr3_end_motif(seq, locus, species):
    #     seq = seq + conserved_aa
    #     if junction_matching_regex.match(seq):
    #         return seq
    #
    # if junction_matching_regex.match(seq):
    #     assert False, "this scenario shouldn happen"
    #     return seq
    #
    # # if any trimming failed, this is likely the CDR3 without conserved aa's; return original sequence
    #
    #
    # # TODO remove *  from motifs
    # # TODO shorter motifs????
    # return orig_seq

