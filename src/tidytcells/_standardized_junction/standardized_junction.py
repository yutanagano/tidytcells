from tidytcells._resources import AMINO_ACIDS
from tidytcells.result import Junction
from tidytcells._utils.alignment import *

C_MISMATCH_AAS = {"W", "S", "R", "G", "Y", "F"}
F_MISMATCH_AAS = {"I", "L", "V", "Y", "S", "C"}
W_MISMATCH_AAS = {"C", "G", "L", "R", "S"}
FW_MISMATCH_AAS = F_MISMATCH_AAS.union(W_MISMATCH_AAS)

MIN_ALIGNMENT_SCORES = {"TR": 2, "IG": 1}

MISMATCH_PENALTY = {"TR": -1.5, "IG": -0.75}


class JunctionStandardizer:
    def __init__(
        self,
        seq: str,
        locus: str,
        j_symbol: str,
        v_symbol: str,
        species: str,
        germline_database: dict,
        allele_tree: dict[str, dict[str, dict]],
        enforce_functional_v: bool = True,
        enforce_functional_j: bool = False,
        allow_c_correction: bool = False,
        allow_fw_correction: bool = False,
        max_v_reconstruction: int = 1,
        max_j_reconstruction: int = 1,
    ) -> None:
        self.orig_seq = seq
        self.corrected_seq = seq
        self.locus = locus
        self.j_symbol = j_symbol
        self.v_symbol = v_symbol
        self.enforce_functional_v = enforce_functional_v
        self.enforce_functional_j = enforce_functional_j
        self.allow_c_correction = allow_c_correction
        self.allow_fw_correction = allow_fw_correction
        self.max_v_reconstruction = max_v_reconstruction
        self.max_j_reconstruction = max_j_reconstruction
        self.species = species
        self.germline_database = germline_database
        self.allele_tree = allele_tree

        self.corrected_first_aa = False
        self.corrected_last_aa = False
        self.correction_j_genes = None
        self.correction_v_genes = None
        self.correction_j_score = None
        self.correction_v_score = None

        self.reasons_invalid = []
        self._resolve_juncton()

        self.result = Junction(
            self.orig_seq,
            self.get_reason_why_invalid(),
            self.corrected_seq,
            self.species,
            self.correction_j_genes,
            self.correction_v_genes,
            self.correction_j_score,
            self.correction_v_score,
            self.germline_database,
        )

    def _resolve_juncton(self) -> None:
        self.j_aa_dict = self.get_aa_dict_from_symbol("J")
        self.v_aa_dict = self.get_aa_dict_from_symbol("V")

        self.j_alignments = self.align_j()
        self.v_alignments = self.align_v()

        if len(self.j_alignments) > 0:
            self.correction_j_score = self.j_alignments[0]["score"]
            self.corrected_seq, self.correction_j_genes = self.correct_seq_j_side(
                self.corrected_seq
            )

        if len(self.v_alignments) > 0:
            self.correction_v_score = self.v_alignments[0]["score"]
            self.corrected_seq, self.correction_v_genes = self.correct_seq_v_side(
                self.corrected_seq
            )

        if len(self.corrected_seq) < 6:
            self.reasons_invalid.append("Junction too short")

    def _get_allele_data(self, allele_symbol: str, gene: str) -> dict | None:
        allele_field_subset = ["FR3-IMGT", "V-REGION"] if gene == "V" else ["J-MOTIF", "J-REGION"]

        gene_name, allele_designation = allele_symbol.split("*")
        allele_data = self.allele_tree.get(gene_name, {}).get(allele_designation)
        allele_data =  {field: allele_data[field] for field in allele_field_subset if field in allele_data}

        if gene == "V":
            v_region_full = allele_data.pop("V-REGION")
            allele_data["V-REGION-SHORT"] = v_region_full

            if "FR3-IMGT" in allele_data:
                try:
                    allele_data["V-REGION-SHORT"] = v_region_full[v_region_full.index(allele_data["FR3-IMGT"]):]
                except ValueError:
                    allele_data.pop("FR3-IMGT")

        return allele_data

    def get_aa_dict_from_symbol(self, gene: str) -> dict[str, dict[str, str]]:
        """
        Given the user-provided allele/gene/subgroup symbol (v_symbol, j_symbol), this function returns
        the sequence data for any of the alleles applicable to the given symbol.
        """

        symbol = self.v_symbol if gene == "V" else self.j_symbol
        symbol = self.locus[0:2] if symbol is None else symbol

        aa_dict = {}

        for sub_symbol in symbol.split(","):
            sub_symbol = sub_symbol.strip()

            # if symbol is an allele, return only info for the given allele
            if "*" in sub_symbol:
                if gene not in sub_symbol:
                    self.reasons_invalid.append(
                        f"not a {gene} gene: {sub_symbol} ({self.species})"
                    )
                    continue

                allele_data = self._get_allele_data(sub_symbol, gene)
                if allele_data is not None and has_sequence_data(allele_data):
                    aa_dict.update({sub_symbol: allele_data})
                    continue
                else:
                    self.reasons_invalid.append(
                        f"no sequence information known for {sub_symbol} ({self.species})"
                    )

            # if sub_symbol is less specific than allele, retrieve any alleles that are valid extensions of the given sub_symbol
            enforce_functional = (
                self.enforce_functional_v if gene == "V" else self.enforce_functional_j
            )
            allele_symbols = get_compatible_symbols(
                sub_symbol,
                self.allele_tree,
                gene,
                self.locus,
                enforce_functional,
            )

            # if all alleles for one gene have the same sequence info, collapse the gene to 1 dict item
            aas_per_allele = {
                ext_symbol: self._get_allele_data(ext_symbol, gene)
                for ext_symbol in allele_symbols
            }
            aas_per_gene = collapse_alleles_with_same_aa(aas_per_allele)
            aa_dict.update(aas_per_gene)

        return aa_dict

    def correct_sequencing_err_j_side(
        self,
        best_alignments_orig: list[dict[str, str | int | float]],
        conserved_aa: str,
    ) -> list[dict[str, str | int | float]]:
        """
        Try correcting sequencing error on the J side (replace last amino acid with conserved_aa; 'F' or 'W').
        If the alignments improve, keep corrected sequence + alignments based on the corrected sequence
        Else, keep original sequence + original alignments
        """

        best_score = (
            -1 if len(best_alignments_orig) == 0 else best_alignments_orig[0]["score"]
        )
        corrected_seq = self.orig_seq[:-1] + conserved_aa

        min_j_score = MIN_ALIGNMENT_SCORES[self.locus[0:2]]
        mismatch_penalty = MISMATCH_PENALTY[self.locus[0:2]]
        corr_best_alignments = align_j_regions(
            corrected_seq,
            self.j_aa_dict,
            min_j_score + 1,
            mismatch_penalty,
            allow_mismatches=False,
        )
        keep_alignments = []

        for alignment in corr_best_alignments:
            # Correct for the "+1" that was added by correcting the anchor residue
            alignment["score"] -= 1

            if alignment["score"] > best_score:
                new_seq_len = alignment["j_offset"] + alignment["j_conserved_idx"] + 1
                if len(corrected_seq) == new_seq_len:
                    keep_alignments.append(alignment)

        if len(keep_alignments) > 0:
            self.corrected_seq = corrected_seq
            self.corrected_last_aa = True
            return keep_alignments

        return best_alignments_orig

    def correct_sequencing_err_v_side(
        self, best_alignments_orig: list[dict[str, str | int | float]]
    ) -> list[dict[str, str | int | float]]:
        """
        Try correcting sequencing error on the V side (replace first amino acid with 'C').
        If the alignments improve, keep corrected sequence + alignments based on the corrected sequence
        Else, keep original sequence + original alignments
        """

        best_score = (
            -1 if len(best_alignments_orig) == 0 else best_alignments_orig[0]["score"]
        )
        corrected_seq = "C" + self.orig_seq[1:]

        min_v_score = MIN_ALIGNMENT_SCORES[self.locus[0:2]]
        mismatch_penalty = MISMATCH_PENALTY[self.locus[0:2]]
        corr_best_alignments = align_v_regions(
            corrected_seq,
            self.v_aa_dict,
            min_v_score + 1,
            mismatch_penalty,
            allow_mismatches=False,
        )
        keep_alignments = []

        for alignment in corr_best_alignments:
            # Correct for the "+1" that was added by correcting the anchor residue
            alignment["score"] -= 1

            if alignment["score"] > best_score:
                # only allow alignments that keep the same length
                if alignment["v_conserved_idx"] == alignment["v_offset"]:
                    keep_alignments.append(alignment)

        if len(keep_alignments) > 0:
            self.corrected_seq = corrected_seq
            self.corrected_first_aa = True
            return keep_alignments

        return best_alignments_orig

    def align_j(self) -> list[dict[str, str | int | float]]:
        """
        Compute alignments for each sequence in self.j_aa_dict, keep only the best alignments
        """
        min_j_score = MIN_ALIGNMENT_SCORES[self.locus[0:2]]
        mismatch_penalty = MISMATCH_PENALTY[self.locus[0:2]]
        best_alignments = align_j_regions(
            self.orig_seq,
            self.j_aa_dict,
            min_j_score,
            mismatch_penalty=mismatch_penalty,
            allow_mismatches=True,
        )

        if self.allow_fw_correction and self.orig_seq[-1] in F_MISMATCH_AAS:
            best_alignments = self.correct_sequencing_err_j_side(
                best_alignments, conserved_aa="F"
            )

        if self.allow_fw_correction and self.orig_seq[-1] in W_MISMATCH_AAS:
            best_alignments = self.correct_sequencing_err_j_side(
                best_alignments, conserved_aa="W"
            )

        if len(best_alignments) == 0:
            err = "J alignment unsuccessful"

            if self.j_symbol is not None:
                err += f" for J symbol {self.j_symbol} ({self.species})"
            else:
                err += f" for any J symbol ({self.species})"

            if len(self.j_aa_dict) == 0:
                err += ": no known sequence information"

            self.reasons_invalid.append(err)

        return best_alignments

    def align_v(self) -> list[dict[str, str | int | float]]:
        """
        Compute alignments for each sequence in self.v_aa_dict, keep only the best alignments
        """
        min_v_score = MIN_ALIGNMENT_SCORES[self.locus[0:2]]
        mismatch_penalty = MISMATCH_PENALTY[self.locus[0:2]]
        best_alignments = align_v_regions(
            self.orig_seq,
            self.v_aa_dict,
            min_v_score,
            mismatch_penalty,
            allow_mismatches=True,
        )

        if self.allow_c_correction and self.orig_seq[0] in C_MISMATCH_AAS:
            best_alignments = self.correct_sequencing_err_v_side(best_alignments)

        if len(best_alignments) == 0:
            err = "V alignment unsuccessful"

            if self.v_symbol is not None:
                err += f" for V symbol {self.v_symbol} ({self.species})"
            else:
                err += " for any V symbol"

            if len(self.v_aa_dict) == 0:
                err += ": no known sequence information"

            self.reasons_invalid.append(err)

        return best_alignments

    def get_match_gene_repr(
        self, seq: str, corrected_seq_with_gene: list[tuple[str, str]]
    ) -> str:
        """
        The gene string representation concatenates all best-scoring V or J genes that match the sequence
        """

        gene_symbols = sorted(
            {
                gene_symbol
                for corrected_seq, gene_symbol in corrected_seq_with_gene
                if corrected_seq == seq
            }
        )

        return ", ".join(gene_symbols)

    def correct_seq_j_side(self, seq: str) -> tuple[str, str | None]:
        """
        Compute corrected sequence for J side

        For each alginment in self.j_alignments, compute the corrected sequence.
        Only return a result if
          1. There exists an alignment matching the sequence perfectly (no correction) or
          2. An unambiguous correction can be determined (no disagreement between alignments)
        """

        corrected_seq_with_j_gene = []

        alignment_too_long = []

        for alignment_details in self.j_alignments:
            j_conserved_idx = alignment_details["j_conserved_idx"]
            j_offset = alignment_details["j_offset"]

            new_seq_len = j_offset + j_conserved_idx + 1

            if len(seq) == new_seq_len:
                corrected_seq_with_j_gene.append((seq, alignment_details["gene"]))

            elif len(seq) > new_seq_len:
                if new_seq_len > 1:
                    corrected_seq_with_j_gene.append(
                        (seq[:new_seq_len], alignment_details["gene"])
                    )

            elif len(seq) < new_seq_len:
                reconstruction_length = new_seq_len - len(seq)

                end_j_idx = j_conserved_idx + 1
                start_j_idx = end_j_idx - reconstruction_length
                reconstructed_aas = alignment_details["j_region"][start_j_idx:end_j_idx]

                if all(aa in AMINO_ACIDS for aa in reconstructed_aas):
                    if reconstruction_length <= self.max_j_reconstruction:
                        corrected_seq_with_j_gene.append(
                            (seq + reconstructed_aas, alignment_details["gene"])
                        )
                    else:
                        alignment_too_long.append(
                            f"{alignment_details['gene']} ({self.species}): Alignment successful but reconstruction too long ({reconstruction_length})"
                        )

        corrected_seqs = {
            corrected_seq for corrected_seq, j_gene in corrected_seq_with_j_gene
        }

        # if the perfect match exists, keep only this
        if seq in corrected_seqs:
            return seq, self.get_match_gene_repr(seq, corrected_seq_with_j_gene)

        if len(corrected_seqs) > 1:
            # When Junction-reconstruction results are ambiguous (different J's)
            # - prefer the junctions which add 1 aa (conserved anchor), if there are any
            # - if the only difference is the last aa, prefer to keep only canonical F/W anchors
            if len(seq) + 1 in {len(s) for s in corrected_seqs}:
                corrected_seqs = {s for s in corrected_seqs if len(s) == len(seq) + 1}

            if len({s[:-1] for s in corrected_seqs}) == 1:
                canonical_seqs = {s for s in corrected_seqs if s[-1] in ("F", "W")}

                if len(canonical_seqs) > 0:
                    corrected_seqs = canonical_seqs

        if len(corrected_seqs) == 0:
            self.reasons_invalid.append(f"J side reconstruction unsuccessful")
            self.reasons_invalid.extend(alignment_too_long)
            return seq, None

        if len(corrected_seqs) > 1:
            self.reasons_invalid.append(
                f"J side reconstruction ambiguous: {corrected_seqs}"
            )
            return seq, None

        corrected_seq = corrected_seqs.pop()
        return corrected_seq, self.get_match_gene_repr(
            corrected_seq, corrected_seq_with_j_gene
        )

    def correct_seq_v_side(self, seq: str) -> tuple[str, str | None]:
        """
        Compute corrected sequence for V side

        For each alginment in self.v_alignments, compute the corrected sequence.
        Only return a result if
          1. There exists an alignment matching the sequence perfectly (no correction) or
          2. An unambiguous correction can be determined (no disagreement between alignments)
        """
        corrected_seq_with_v_gene = []
        alignment_too_long = []

        for alignment_details in self.v_alignments:
            v_conserved_idx = alignment_details["v_conserved_idx"]
            v_offset = alignment_details["v_offset"]

            if v_conserved_idx == v_offset:
                corrected_seq_with_v_gene.append((seq, alignment_details["gene"]))

            if v_conserved_idx > v_offset:
                n_aas_to_trim = v_conserved_idx - v_offset
                new_seq = seq[n_aas_to_trim:]
                corrected_seq_with_v_gene.append((new_seq, alignment_details["gene"]))

            if v_conserved_idx < v_offset:
                reconstructed_aas = alignment_details["v_region_short"][v_conserved_idx:][
                    :v_offset
                ]

                if all(aa in AMINO_ACIDS for aa in reconstructed_aas):
                    reconstruction_length = len(reconstructed_aas)

                    if reconstruction_length <= self.max_v_reconstruction:
                        new_seq = reconstructed_aas + seq
                        corrected_seq_with_v_gene.append(
                            (new_seq, alignment_details["gene"])
                        )
                    else:
                        alignment_too_long.append(
                            f"{alignment_details['gene']} ({self.species}): Alignment successful but reconstruction too long ({reconstruction_length})"
                        )

        corrected_seq_with_v_gene = [
            (corrected_seq, v_gene)
            for corrected_seq, v_gene in corrected_seq_with_v_gene
            if corrected_seq.startswith("C")
        ]
        corrected_seqs = {
            corrected_seq for corrected_seq, v_gene in corrected_seq_with_v_gene
        }

        # if the perfect match exists, keep only this
        if seq in corrected_seqs:
            return seq, self.get_match_gene_repr(seq, corrected_seq_with_v_gene)

        # prefer the junction which only adds the conserved C, if there is one
        if "C" + seq in corrected_seqs:
            return "C" + seq, self.get_match_gene_repr(
                "C" + seq, corrected_seq_with_v_gene
            )

        if len(corrected_seqs) == 0:
            self.reasons_invalid.append(f"V side reconstruction unsuccessful")
            self.reasons_invalid.extend(alignment_too_long)
            return seq, None

        if len(corrected_seqs) > 1:
            self.reasons_invalid.append(
                f"V side reconstruction ambiguous: {corrected_seqs}"
            )
            return seq, None

        corrected_seq = corrected_seqs.pop()
        return corrected_seq, self.get_match_gene_repr(
            corrected_seq, corrected_seq_with_v_gene
        )

    def get_reason_why_invalid(
        self,
    ) -> str | None:  # optional: enforce_complete / no reconstruction etc
        """
        If the CDR3 cannot be standardized (it is invalid), this method returns a string outlining the reason why (incomplete on the left side, right side, etc).
        Returns None if standardization was successful.
        """
        if len(self.reasons_invalid) == 0:
            return None

        return "; ".join(self.reasons_invalid)