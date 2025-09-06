import logging
import re

from tidytcells._resources import HOMOSAPIENS_TR_AA_SEQUENCES, HOMOSAPIENS_IG_AA_SEQUENCES, MUSMUSCULUS_TR_AA_SEQUENCES


logger = logging.getLogger(__name__)


def get_conserved_v_motifs_from_aa_dict(aa_dict, locus):
    motifs = set()

    locus = locus if locus is not None else ""

    for allele, seq_dict in aa_dict.items():
        if allele.startswith(locus):
            if "FR3-IMGT" in seq_dict:
                if seq_dict["FR3-IMGT"].endswith("C"):
                    motifs.add(seq_dict["FR3-IMGT"][-4:])

    return motifs

def get_conserved_v_motifs(locus, species):
    # Return all motifs of the last 4 amino acids in a V gene, ending at the conserved C start of the CDR3
    v_motifs = set()

    if species == "homosapiens":
        if locus is None:
            v_motifs = set()
            v_motifs.update(get_conserved_v_motifs_from_aa_dict(HOMOSAPIENS_TR_AA_SEQUENCES, locus))
            v_motifs.update(get_conserved_v_motifs_from_aa_dict(HOMOSAPIENS_IG_AA_SEQUENCES, locus))
        elif locus.startswith("TR"):
            v_motifs = get_conserved_v_motifs_from_aa_dict(HOMOSAPIENS_TR_AA_SEQUENCES, locus)
        elif locus.startswith("IG"):
            v_motifs = get_conserved_v_motifs_from_aa_dict(HOMOSAPIENS_IG_AA_SEQUENCES, locus)
    elif species == "musmusculus":
        if locus is not None and locus.startswith("IG"):
            raise ValueError(f"Trimming is not supported for locus {locus} and species {species}.")
        v_motifs = get_conserved_v_motifs_from_aa_dict(MUSMUSCULUS_TR_AA_SEQUENCES, locus)
    else:
        raise ValueError(f"Trimming is not supported for species {species}.")

    return v_motifs


def look_for_v_motif_start_recursive(orig_seq, v_motifs):
    v_motifs = {motif[1:] for motif in v_motifs}

    if v_motifs == {"C"} or v_motifs == {""} or v_motifs == set():
        return orig_seq

    for motif in v_motifs:
        if orig_seq.startswith(motif):
            return orig_seq.replace(motif, "C", 1)

    return look_for_v_motif_start_recursive(orig_seq, v_motifs)


def trim_junction_start(orig_seq, locus, species):
    v_motifs = get_conserved_v_motifs(locus, species)

    # case 1: long sequence provided (>4 extra aa's at start)
    #   look for conserved V motif of length 4 *anywhere* in the sequence. if found, trim
    #   due to the length, this rarely occurs randomly in CDR3
    for motif in v_motifs:
        if motif in orig_seq:
            new_seq = "C" + orig_seq.split(motif)[1]
            if len(new_seq) >= 6:
                return new_seq

    # case 2: few extra amino acids at start, specific known motifs
    #   trim these *only* if they appear at the start
    return look_for_v_motif_start_recursive(orig_seq, v_motifs)


def trim_junction(orig_seq, junction_matching_regex, locus, species):
    # do not trim if it is already a valid CDR3 junction
    if junction_matching_regex.match(orig_seq):
        return orig_seq

    seq = trim_junction_start(orig_seq, locus, species)

    junction_matching_regex_open_ended = re.compile(junction_matching_regex.pattern[:-1])
    match_obj = junction_matching_regex_open_ended.match(seq)

    if match_obj:
        seq = match_obj.group(0)

        if len(seq) >= 6:
            return seq

    # if no match (= no trailing conserved aa), or the resulting sequence is too short:
    # this is more likely to be an IMGT-CDR3; return untrimmed sequence so leading/trailing aa's can be added
    return orig_seq


