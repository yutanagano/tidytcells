import logging

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

    assert motif_type in ("V-MOTIF", "J-MOTIF", "V-CDR3-START", "J-CDR3-END"), f"Unknown motif type: {motif_type}"

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


def startswith_cdr3_motif(seq, locus, species, include_conserved, max_motif_len=2):
    assert type(include_conserved) == bool

    start_motifs = get_conserved_motifs(locus, species, "V-CDR3-START")
    start_motifs.discard("C")

    if max_motif_len is not None:
        start_motifs = {motif[:max_motif_len] if len(motif) > max_motif_len else motif for motif in start_motifs}

    if not include_conserved:
        start_motifs = {motif[1:] for motif in start_motifs}

    for motif in start_motifs:
        if seq.startswith(motif):
            return True

    return False # todo one of the motifs is literally 'C'.. defeats the purpose?


def endswith_cdr3_motif(seq, locus, species, include_conserved=False, max_motif_len=2):
    assert type(include_conserved) == bool

    end_motifs = get_conserved_motifs(locus, species, "J-CDR3-END")

    if max_motif_len is not None:
        end_motifs = {motif[-max_motif_len:] if len(motif) > max_motif_len else motif for motif in end_motifs}

    end_motifs.discard("")

    if not include_conserved:
        end_motifs = {motif[:-1] for motif in end_motifs}

    for motif in end_motifs:
        if seq.endswith(motif):
            return True

    return False


def minimal_trim_junction_start(orig_seq, locus, species):
    # if junction / or CDR3
    if startswith_cdr3_motif(orig_seq, locus, species, include_conserved=True, max_motif_len=None) or startswith_cdr3_motif(orig_seq, locus, species, include_conserved=False, max_motif_len=None):
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
    if endswith_cdr3_motif(orig_seq, locus, species, include_conserved=True, max_motif_len=None) or endswith_cdr3_motif(orig_seq, locus, species, include_conserved=False, max_motif_len=None):
        return orig_seq

    # if orig_seq[-1] in conserved_aas or endswith_cdr3_end_motif(orig_seq, locus, species):
    #     return orig_seq

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


def is_valid_junction(seq, junction_matching_regex, locus, species, check_motifs=True):
    if junction_matching_regex.match(seq):

        if check_motifs:
            if startswith_cdr3_motif(seq, locus, species, include_conserved=True) and endswith_cdr3_motif(seq, locus, species, include_conserved=True):
                return True
        else:
            return True

    return False


def is_valid_cdr3(seq, locus, species, min_len=4):
    if len(seq) >= min_len:
        if startswith_cdr3_motif(seq, locus, species, include_conserved=False) and endswith_cdr3_motif(seq, locus, species, include_conserved=False):
            return True

    return False


def _get_conserved_aa(gene_seqs_dict):
    if "J-PHE" in gene_seqs_dict and gene_seqs_dict["J-PHE"] == "F":
        return "F"
    elif "J-TRP" in gene_seqs_dict and gene_seqs_dict["J-TRP"] == "W":
        return "W"
    elif "J-CYS" in gene_seqs_dict and gene_seqs_dict["J-CYS"] == "C":
        return "C"
    elif "J-VAL" in gene_seqs_dict and gene_seqs_dict["J-VAL"] == "V":
        return "V"



def process_junction(orig_seq, junction_matching_regex, conserved_aa, locus, species, trimming=True, check_motifs=True):
    # cases are checked in order of most likely occurrence + minimal interventions
    # case 1: is junction with leading/trailing, CAAAAF (already valid)
    # case 2: is CDR3 without leading/trailing, AAAA (reconstruct leading/trailing)
    # case 3: is junction with excess on 1 or both sides (fully recoverable substring) xxxCAAAAFxxx / CAAAAFxxx / xxxCAAAAF
    # case 4: is CDR3 with excess on 1 or both sides (mostly substring, reconstruct leading/trailing) AAAAFxxx / xxxCAAAA
    # case 5: is partial CDR3 / does not contain valid motifs (needs reconstruction?? only for TCR)

    # case 1: junction
    if is_valid_junction(orig_seq, junction_matching_regex, locus, species, check_motifs):
        # print("valid junction:", orig_seq) # todo for the cdr3 motifs keep with conserved aa
        return orig_seq

    # case 2: cdr3
    if check_motifs and is_valid_cdr3(orig_seq, locus, species):
        # todo retrieve precise conserved aa if all_conserved_aas is ambiguous?, this one is just preferred
        # print("valid cdr3 to junction:", orig_seq, ":", "C" + orig_seq + conserved_aa)

        return "C" + orig_seq + conserved_aa

    all_conserved_aas = get_expected_conserved_aas(junction_matching_regex)  # todo if all_conserved_aas >1, determine based on J motif

    # trimming
    if trimming:
        seq = minimal_trim_junction_start(orig_seq, locus, species)
        seq = minimal_trim_junction_end(seq, locus, species, all_conserved_aas)
    else:
        seq = orig_seq

    # case 3
    if is_valid_junction(seq, junction_matching_regex, locus, species, check_motifs):
        # print("valid junction after trimming:", orig_seq, ":", seq)
        return seq

    # case 4 (only add conserved aas to sides that were not trimmed)
    start_was_trimmed = not orig_seq.startswith(seq)
    end_was_trimmed = not orig_seq.endswith(seq)

    # only 'add' to sites that were not trimmed
    #   if the current start/end is not a valid conserved motif
    #     but adding a C / FW makes it a conserved motif
    #       then add C / FW

    if not start_was_trimmed:
        if not startswith_cdr3_motif(seq, locus, species, include_conserved=True):
            if startswith_cdr3_motif("C" + seq, locus, species, include_conserved=True):
                seq = "C" + seq

    if not end_was_trimmed:
        if not endswith_cdr3_motif(seq, locus, species, include_conserved=True):
            if endswith_cdr3_motif(seq + conserved_aa, locus, species, include_conserved=True):
                seq = seq + conserved_aa
            else: # special case: check other possible conserved end aa's
                for aa in all_conserved_aas:
                    if aa != conserved_aa and endswith_cdr3_motif(seq + aa, locus, species, include_conserved=True):
                        seq = seq + aa
                        break


    # if startswith_cdr3_motif(seq, locus, species, include_conserved=False) and not start_was_trimmed:
    #     seq = "C" + seq

    # if endswith_cdr3_motif(seq, locus, species, include_conserved=False) and not end_was_trimmed:
    #     seq = seq + conserved_aa

    if is_valid_junction(seq, junction_matching_regex, locus, species, check_motifs):
        # print("valid junction after trimming and restoring:", orig_seq, ":", seq)
        return seq

    # assert False, "incomplete or invalid CDR3?"
    # case 5: invalid CDR3

    if startswith_cdr3_motif(seq, locus, species, include_conserved=True) and endswith_cdr3_motif(seq, locus, species, include_conserved=True):
        pass

    # print(locus, orig_seq, ":", seq, f"start={startswith_cdr3_motif(seq, locus, species, include_conserved=True)}", f"end={endswith_cdr3_motif(seq, locus, species, include_conserved=True)}", sep="\t")

    return orig_seq

   # 'AVKDARLMFG' -> 'CAVKDARLMF'
   # 'ARHRNWLFDY' / 'CARHRNWLFDY'
   # 'AVLNTGGFKTI' / 'CAVLNTGGFKTI' -> NTGGFKTI
   # 'AAFDDKIIFG' / 'CAAFDDKIIF'
   # 'AANNARLMFG' / 'CAANNARLMF'
   # 'AASASKLIFG' / 'CAASASKLIF'
   # 'AASGTGAGSYQLT' / 'CAASGTGAGSYQLT'
   # 'AASSLYGQNFV' / 'CAASSLYGQNFV'


    #
    # print("incomplete or invalid CDR3:", orig_seq, f": {seq})" if seq != orig_seq else "")
    # return None
    #
    #
    #
    #
    #


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


