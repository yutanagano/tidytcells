from abc import ABC, abstractmethod
from typing import Optional, Dict

from tidytcells._standardized_results.StandardizedResult import StandardizedJunctionResult
from tidytcells._utils.alignment import *


C_MISMATCH_AAS = {"W", "S", "R", "G", "Y", "F"}
F_MISMATCH_AAS = {"I", "L", "V", "Y", "S", "C"}
W_MISMATCH_AAS = {"C", "G", "L", "R", "S"}
FW_MISMATCH_AAS = F_MISMATCH_AAS.union(W_MISMATCH_AAS)

class StandardizedJunction(ABC):
    """
    Abstract base standardizer class.
    """

    @property
    @abstractmethod
    def _sequence_dictionary(self) -> Dict[str, Dict]:
        pass

    def __init__(self, seq: str, locus: str, j_symbol: str, v_symbol: str,
                 enforce_functional_v: bool = True, enforce_functional_j: bool = False,
                 allow_c_correction: bool = False, allow_fw_correction: bool = False,
                 allow_v_reconstruction: bool = False, allow_j_reconstruction: bool = False,
                 mismatch_penalty: float = -1.5, max_v_mismatches: int = 0, max_j_mismatches: int = 1,
                 min_j_score: int = 2, min_v_score: int = 1) -> None:
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
        self.mismatch_penalty = mismatch_penalty
        self.max_v_mismatches = max_v_mismatches
        self.max_j_mismatches = max_j_mismatches
        self.min_j_score = min_j_score
        self.min_v_score = min_v_score
        self.corrected_first_aa = False
        self.corrected_last_aa = False

        self.reasons_invalid = []
        self._resolve_juncton()

        self.result = StandardizedJunctionResult(self.orig_seq, "; ".join(self.reasons_invalid), self.corrected_seq)


    def _resolve_juncton(self):
        self.j_aa_dict = self.get_aa_dict_from_symbol("J")
        self.v_aa_dict = self.get_aa_dict_from_symbol("V")

        self.j_alignments = self.align_j()
        self.v_alignments = self.align_v()

        self.corrected_seq = self.correct_seq_j_side(self.corrected_seq)
        self.corrected_seq = self.correct_seq_v_side(self.corrected_seq)

        if len(self.corrected_seq) < 4:
            self.reasons_invalid.append("junction too short")

    def get_aa_dict_from_symbol(self, gene) -> dict:
        symbol = self.locus[0:2]

        if gene == "J" and self.j_symbol is not None:
            symbol = self.j_symbol

        if gene == "V" and self.v_symbol is not None:
            symbol = self.v_symbol

        if "*" in symbol:
            if gene not in symbol:
                self.reasons_invalid.append(f"not a {gene} gene: {symbol}")
                return dict()

            if symbol in  self._sequence_dictionary:
                return {symbol: self._sequence_dictionary[symbol]}
            else:
                self.reasons_invalid.append("no sequence information known for " + symbol)

        enforce_functional = self.enforce_functional_v if gene == "V" else self.enforce_functional_j

        allele_symbols = get_compatible_symbols(symbol, self._sequence_dictionary, gene, self.locus, enforce_functional)
        aas_per_allele = {ext_symbol: self._sequence_dictionary[ext_symbol] for ext_symbol in allele_symbols}
        aas_per_gene = collapse_aa_dict_per_gene(aas_per_allele)

        return aas_per_gene

    def align_j(self):
        best_alignments = align_j_regions(self.orig_seq, self.j_aa_dict, self.min_j_score, self.mismatch_penalty, self.max_j_mismatches)

        if self.allow_fw_correction and self.orig_seq[-1] in FW_MISMATCH_AAS:

            # todo fill in part with F/W
            #   one solution:

            # TRY F FIRST # todo add some locus restriction??
            best_score = -1 if len(best_alignments) == 0 else best_alignments[0]["score"]
            corrected_seq = self.orig_seq[:-1] + "F"
            corr_best_alignments = align_j_regions(corrected_seq, self.j_aa_dict, self.min_j_score+1,
                                                   self.mismatch_penalty, max_mismatches=0)
            corr_best_score = -1 if len(corr_best_alignments) == 0 else corr_best_alignments[0]["score"]

            if corr_best_score > best_score:
                self.corrected_seq = corrected_seq
                self.corrected_last_aa = True
                best_alignments = corr_best_alignments

            # THEN TRY W (basically, if F doesnt align) # todo add some locus restriction
            best_score = -1 if len(best_alignments) == 0 else best_alignments[0]["score"]
            corrected_seq = self.orig_seq[:-1] + "W"
            corr_best_alignments = align_j_regions(corrected_seq, self.j_aa_dict, self.min_j_score+1,
                                                   self.mismatch_penalty, max_mismatches=0)
            corr_best_score = -1 if len(corr_best_alignments) == 0 else corr_best_alignments[0]["score"]

            if corr_best_score > best_score:
                self.corrected_seq = corrected_seq
                self.corrected_last_aa = True
                best_alignments = corr_best_alignments # todo check all these cases

        if len(best_alignments) == 0:
            err = "J alignment unsuccessful"
            if self.j_symbol is not None:
                err += " for J symbol " + self.j_symbol

            if len(self.j_aa_dict) == 0:
                err += ": no known sequence information"

            self.reasons_invalid.append(err)

        return best_alignments

    def align_v(self):
        best_alignments = align_v_regions(self.orig_seq, self.v_aa_dict, self.min_v_score, self.mismatch_penalty, self.max_v_mismatches)

        if self.allow_c_correction and self.orig_seq[0] in C_MISMATCH_AAS:
            best_score = -1 if len(best_alignments) == 0 else best_alignments[0]["score"]
            corrected_seq = "C" + self.orig_seq[1:]

            corr_best_alignments = align_v_regions(corrected_seq, self.v_aa_dict, self.min_v_score+1,
                                                   self.mismatch_penalty, 0)
            corr_best_score = -1 if len(corr_best_alignments) == 0 else corr_best_alignments[0]["score"]

            if corr_best_score > best_score:
                self.corrected_seq = corrected_seq
                self.corrected_first_aa = True
                best_alignments = corr_best_alignments

        if len(best_alignments) == 0:
            err = "V alignment unsuccessful"
            if self.v_symbol is not None:
                err += " for V symbol " + self.v_symbol

            if len(self.v_aa_dict) == 0:
                err += ": no known sequence information"

            self.reasons_invalid.append(err)

        return best_alignments

    def correct_seq_j_side(self, seq):
        results = set()

        for alignment_details in self.j_alignments:
            j_conserved_idx = alignment_details["j_conserved_idx"]
            j_offset = alignment_details["j_offset"]

            new_seq_len = j_offset + j_conserved_idx + 1

            # If any alignment matches perfectly, don't look any further
            if len(seq) == new_seq_len:
                return seq

            elif len(seq) > new_seq_len:
                results.add(seq[:new_seq_len])

            elif len(seq) < new_seq_len:
                reconstruction_length = new_seq_len - len(seq)

                if reconstruction_length <= 1 or self.allow_j_reconstruction:
                    end_j_idx = j_conserved_idx + 1
                    start_j_idx = end_j_idx - reconstruction_length

                    reconstructed_aas = alignment_details["j_region"][start_j_idx:end_j_idx]
                    new_seq = seq + reconstructed_aas
                    results.add(new_seq)

        if len(results) > 1:
            # When Junction-reconstruction results are ambiguous (different J's)
            # - prefer the junctions which add 1 aa (conserved anchor), if there are any
            # - without J symbol, prefer to keep only canonical F/W anchors

            if len(seq) + 1 in {len(s) for s in results}:
                results = {s for s in results if len(s) == len(seq) + 1}

            if self.j_symbol is None:
                if len(results) > 1:
                    results = {s for s in results if s[-1] in ("F", "W")}

            if len(results) > 1:
                self.reasons_invalid.append(f"J side reconstruction ambiguous: {results}")
                return seq

        if len(results) == 0:
            self.reasons_invalid.append(f"J side reconstruction unsuccessful.")
            return seq

        return results.pop()

    def correct_seq_v_side(self, seq):
        results = set()

        for alignment_details in self.v_alignments:
            v_conserved_idx = alignment_details["v_conserved_idx"]
            v_offset = alignment_details["v_offset"]

            if v_conserved_idx == v_offset:
                return seq

            if v_conserved_idx > v_offset:
                n_aas_to_trim = v_conserved_idx - v_offset
                new_seq = seq[n_aas_to_trim:]
                results.add(new_seq)

            if v_conserved_idx < v_offset:
                reconstructed_aas = alignment_details["v_region"][v_conserved_idx:][:v_offset]

                if len(reconstructed_aas) <= 1 or self.allow_v_reconstruction:
                    new_seq = reconstructed_aas + seq
                    results.add(new_seq)

        results = {seq for seq in results if seq.startswith("C")}

        if "C" + seq in results:
            return "C" + seq

        if len(results) == 0:
            self.reasons_invalid.append(f"V side reconstruction unsuccessful.")
            return seq

        if len(results) > 1:
            self.reasons_invalid.append(f"V side reconstruction ambiguous: {results}")
            return seq

        return results.pop()

    def get_reason_why_invalid(self) -> Optional[str]: # optional: enforce_complete / no reconstruction etc
        """
        If the CDR3 cannot be standardized (it is invalid), this method returns a string outlining the reason why (incomplete on the left side, right side, etc).
        Returns None if standardisation was successful.
        """
        # todo optional: check if the same locus is used for V and J.
        #  however, the difficulty is that multiple alignments (of different loci) are possible

        if len(self.reasons_invalid) == 0:
            return None

        return ", ".join(self.reasons_invalid)

    # def compile(self, region: str = "junction") -> str:
    #     """
    #     Compile a complete string representation of the gene.
    #     The argument given to precision will determine the amount of specificity given in the compiled string.
    #     """
    #     if self.get_reason_why_invalid() is None:
    #         if region.lower() == "junction":
    #             return self.corrected_seq
    #
    #         if region.lower() == "cdr3":
    #             return self.corrected_seq[1:-1]
    #
