from abc import ABC, abstractmethod
from typing import Dict, Optional

from tidytcells._utils.result import Junction
from tidytcells._utils.alignment import *


C_MISMATCH_AAS = {"W", "S", "R", "G", "Y", "F"}
F_MISMATCH_AAS = {"I", "L", "V", "Y", "S", "C"}
W_MISMATCH_AAS = {"C", "G", "L", "R", "S"}
FW_MISMATCH_AAS = F_MISMATCH_AAS.union(W_MISMATCH_AAS)

MIN_J_SCORE = 1
MIN_V_SCORE = 1
MAX_J_MISMATCHES = 1
MAX_V_MISMATCHES = 0
MISMATCH_PENALTY = -1.5


class JunctionStandardizer(ABC):
    """
    Abstract base standardizer class.
    """

    @property
    @abstractmethod
    def _species(self) -> str:
        pass

    @property
    @abstractmethod
    def _sequence_dictionary(self) -> Dict[str, Dict]:
        pass

    def __init__(self, seq: str, locus: str, j_symbol: str, v_symbol: str,
                 enforce_functional_v: bool = True, enforce_functional_j: bool = False,
                 allow_c_correction: bool = False, allow_fw_correction: bool = False,
                 allow_v_reconstruction: bool = False, allow_j_reconstruction: bool = False) -> None:
        self.orig_seq = seq
        self.corrected_seq = seq
        self.locus = locus
        self.j_symbol = j_symbol
        self.v_symbol = v_symbol
        self.enforce_functional_v = enforce_functional_v
        self.enforce_functional_j = enforce_functional_j
        self.allow_c_correction = allow_c_correction
        self.allow_fw_correction = allow_fw_correction
        self.allow_v_reconstruction = allow_v_reconstruction
        self.allow_j_reconstruction = allow_j_reconstruction
        self.corrected_first_aa = False
        self.corrected_last_aa = False

        self.reasons_invalid = []
        self._resolve_juncton()

        self.result = Junction(self.orig_seq, self.get_reason_why_invalid(), self.corrected_seq, self._species)


    def _resolve_juncton(self):
        self.j_aa_dict = self.get_aa_dict_from_symbol("J")
        self.v_aa_dict = self.get_aa_dict_from_symbol("V")

        self.j_alignments = self.align_j()
        self.v_alignments = self.align_v()

        self.corrected_seq = self.correct_seq_j_side(self.corrected_seq)
        self.corrected_seq = self.correct_seq_v_side(self.corrected_seq)

        if len(self.corrected_seq) < 6:
            self.reasons_invalid.append("junction too short")

    def get_aa_dict_from_symbol(self, gene) -> dict:
        '''
        Given the user-provided allele/gene/subgroup symbol (v_symbol, j_symbol), this function returns
        the sequence_dictionary information for any of the alleles applicable to the given symbol.
        '''

        symbol = self.locus[0:2]

        if gene == "J" and self.j_symbol is not None:
            symbol = self.j_symbol

        if gene == "V" and self.v_symbol is not None:
            symbol = self.v_symbol

        # if symbol is an allele, return only info for the given allele
        if "*" in symbol:
            if gene not in symbol:
                self.reasons_invalid.append(f"not a {gene} gene: {symbol}")
                return dict()

            if symbol in  self._sequence_dictionary:
                return {symbol: self._sequence_dictionary[symbol]}
            else:
                self.reasons_invalid.append("no sequence information known for " + symbol)

        # if symbol is less specific than allele, retrieve any alleles that are valid extensions of the given symbol
        enforce_functional = self.enforce_functional_v if gene == "V" else self.enforce_functional_j
        allele_symbols = get_compatible_symbols(symbol, self._sequence_dictionary, gene, self.locus, enforce_functional)

        # if all alleles for one gene have the same sequence info, collapse the gene to 1 dict item
        aas_per_allele = {ext_symbol: self._sequence_dictionary[ext_symbol] for ext_symbol in allele_symbols}
        aas_per_gene = collapse_aa_dict_per_gene(aas_per_allele)

        return aas_per_gene

    def correct_sequencing_err_j_side(self, best_alignments_orig, conserved_aa):
        '''
        Try correcting sequencing error on the J side (replace last amino acid with conserved_aa; 'F' or 'W').
        If the alignments improve, keep corrected sequence + alignments based on the corrected sequence
        Else, keep original sequence + original alignments
        '''

        best_score = -1 if len(best_alignments_orig) == 0 else best_alignments_orig[0]["score"]
        corrected_seq = self.orig_seq[:-1] + conserved_aa

        corr_best_alignments = align_j_regions(corrected_seq, self.j_aa_dict, MIN_J_SCORE + 1,
                                               MISMATCH_PENALTY, max_mismatches=0)
        keep_alignments = []

        for alignment in corr_best_alignments:
            if alignment["score"] > best_score:
                new_seq_len = alignment["j_offset"] + alignment["j_conserved_idx"] + 1
                if len(corrected_seq) == new_seq_len:
                    keep_alignments.append(alignment)

        if len(keep_alignments) > 0:
            self.corrected_seq = corrected_seq
            self.corrected_last_aa = True
            return keep_alignments

        return best_alignments_orig


    def correct_sequencing_err_v_side(self, best_alignments_orig):
        '''
        Try correcting sequencing error on the V side (replace first amino acid with 'C').
        If the alignments improve, keep corrected sequence + alignments based on the corrected sequence
        Else, keep original sequence + original alignments
        '''

        best_score = -1 if len(best_alignments_orig) == 0 else best_alignments_orig[0]["score"]
        corrected_seq = "C" + self.orig_seq[1:]

        corr_best_alignments = align_v_regions(corrected_seq, self.v_aa_dict, MIN_V_SCORE + 1,
                                               MISMATCH_PENALTY, 0)
        keep_alignments = []

        for alignment in corr_best_alignments:
            if alignment["score"] > best_score:
                # only allow alignments that keep the same length
                if alignment["v_conserved_idx"] == alignment["v_offset"]:
                    keep_alignments.append(alignment)

        if len(keep_alignments) > 0:
            self.corrected_seq = corrected_seq
            self.corrected_first_aa = True
            return keep_alignments

        return best_alignments_orig


    def align_j(self):
        '''
        Compute alignments for each sequence in self.j_aa_dict, keep only the best alignments
        '''

        best_alignments = align_j_regions(self.orig_seq, self.j_aa_dict, MIN_J_SCORE, MISMATCH_PENALTY, MAX_J_MISMATCHES)

        if self.allow_fw_correction and self.orig_seq[-1] in F_MISMATCH_AAS:
            best_alignments = self.correct_sequencing_err_j_side(best_alignments, conserved_aa="F")

        if self.allow_fw_correction and self.orig_seq[-1] in W_MISMATCH_AAS:
            best_alignments = self.correct_sequencing_err_j_side(best_alignments, conserved_aa="W")

        if len(best_alignments) == 0:
            err = "J alignment unsuccessful"
            if self.j_symbol is not None:
                err += " for J symbol " + self.j_symbol

            if len(self.j_aa_dict) == 0:
                err += ": no known sequence information"

            self.reasons_invalid.append(err)

        return best_alignments

    def align_v(self):
        '''
        Compute alignments for each sequence in self.v_aa_dict, keep only the best alignments
        '''

        best_alignments = align_v_regions(self.orig_seq, self.v_aa_dict, MIN_V_SCORE, MISMATCH_PENALTY, MAX_V_MISMATCHES)

        if self.allow_c_correction and self.orig_seq[0] in C_MISMATCH_AAS:
            best_alignments = self.correct_sequencing_err_v_side(best_alignments)

        if len(best_alignments) == 0:
            err = "V alignment unsuccessful"
            if self.v_symbol is not None:
                err += " for V symbol " + self.v_symbol

            if len(self.v_aa_dict) == 0:
                err += ": no known sequence information"

            self.reasons_invalid.append(err)

        return best_alignments

    def correct_seq_j_side(self, seq):
        '''
        Compute corrected sequence for J side

        For each alginment in self.j_alignments, compute the corrected sequence.
        Only return a result if
          1. There exists an alignment matching the sequence perfectly (no correction) or
          2. An unambiguous correction can be determined (no disagreement between alignments)
        '''

        corrected_seqs = set()

        for alignment_details in self.j_alignments:
            j_conserved_idx = alignment_details["j_conserved_idx"]
            j_offset = alignment_details["j_offset"]

            new_seq_len = j_offset + j_conserved_idx + 1

            # If any alignment matches perfectly, don't look any further
            if len(seq) == new_seq_len:
                return seq

            elif len(seq) > new_seq_len:
                corrected_seqs.add(seq[:new_seq_len])

            elif len(seq) < new_seq_len:
                reconstruction_length = new_seq_len - len(seq)

                if reconstruction_length <= 1 or self.allow_j_reconstruction:
                    end_j_idx = j_conserved_idx + 1
                    start_j_idx = end_j_idx - reconstruction_length

                    reconstructed_aas = alignment_details["j_region"][start_j_idx:end_j_idx]
                    corrected_seqs.add(seq + reconstructed_aas)

        if len(corrected_seqs) > 1:
            # When Junction-reconstruction results are ambiguous (different J's)
            # - prefer the junctions which add 1 aa (conserved anchor), if there are any
            # - without J symbol, prefer to keep only canonical F/W anchors

            if len(seq) + 1 in {len(s) for s in corrected_seqs}:
                corrected_seqs = {s for s in corrected_seqs if len(s) == len(seq) + 1}

            if self.j_symbol is None:
                if len(corrected_seqs) > 1:
                    corrected_seqs = {s for s in corrected_seqs if s[-1] in ("F", "W")}

            if len(corrected_seqs) > 1:
                self.reasons_invalid.append(f"J side reconstruction ambiguous: {corrected_seqs}")
                return seq

        if len(corrected_seqs) == 0:
            self.reasons_invalid.append(f"J side reconstruction unsuccessful.")
            return seq

        return corrected_seqs.pop()

    def correct_seq_v_side(self, seq):
        '''
        Compute corrected sequence for V side

        For each alginment in self.v_alignments, compute the corrected sequence.
        Only return a result if
          1. There exists an alignment matching the sequence perfectly (no correction) or
          2. An unambiguous correction can be determined (no disagreement between alignments)
        '''
        corrected_seqs = set()

        for alignment_details in self.v_alignments:
            v_conserved_idx = alignment_details["v_conserved_idx"]
            v_offset = alignment_details["v_offset"]

            if v_conserved_idx == v_offset:
                return seq

            if v_conserved_idx > v_offset:
                n_aas_to_trim = v_conserved_idx - v_offset
                new_seq = seq[n_aas_to_trim:]
                corrected_seqs.add(new_seq)

            if v_conserved_idx < v_offset:
                reconstructed_aas = alignment_details["v_region"][v_conserved_idx:][:v_offset]

                if len(reconstructed_aas) <= 1 or self.allow_v_reconstruction:
                    new_seq = reconstructed_aas + seq
                    corrected_seqs.add(new_seq)

        corrected_seqs = {seq for seq in corrected_seqs if seq.startswith("C")}

        if "C" + seq in corrected_seqs:
            return "C" + seq

        if len(corrected_seqs) == 0:
            self.reasons_invalid.append(f"V side reconstruction unsuccessful.")
            return seq

        if len(corrected_seqs) > 1:
            self.reasons_invalid.append(f"V side reconstruction ambiguous: {corrected_seqs}")
            return seq

        return corrected_seqs.pop()

    def get_reason_why_invalid(self) -> Optional[str]: # optional: enforce_complete / no reconstruction etc
        """
        If the CDR3 cannot be standardized (it is invalid), this method returns a string outlining the reason why (incomplete on the left side, right side, etc).
        Returns None if standardisation was successful.
        """
        if len(self.reasons_invalid) == 0:
            return None

        return "; ".join(self.reasons_invalid)

