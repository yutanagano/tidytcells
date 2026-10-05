import collections
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Tuple
import script_utility


PATH_TO_DATA_DIR = Path("data") / "OGRDB"

LOCI = ("IGH", "IGK", "IGL")

V_REGIONS = (
    ("fwr1", "FR1-IMGT"),
    ("cdr1", "CDR1-IMGT"),
    ("fwr2", "FR2-IMGT"),
    ("cdr2", "CDR2-IMGT"),
    ("fwr3", "FR3-IMGT"),
)

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def main() -> None:
    print("Fetching IG germline sets from OGRDB...")
    paths = tuple(get_germline_set_path(locus) for locus in LOCI)

    print("Converting IG gene sequence data...")
    sequence_data = get_ig_aa_sequence_data(paths)
    script_utility.save_as_json(sequence_data, "homosapiens_ig_aa_sequences_ogrdb.json")


def get_germline_set_path(locus: str) -> Path:
    path = PATH_TO_DATA_DIR / f"Homo_sapiens_{locus}.json"

    if not path.exists():
        download_germline_set(locus, path)

    return path


def download_germline_set(locus: str, path: Path) -> None:
    if shutil.which("download_germline_set") is None:
        raise RuntimeError("download_germline_set not found, install with 'pip install receptor-utils'")

    command = ["download_germline_set", "Homo sapiens", locus, "-f", "AIRRC-JSON"]

    if locus == "IGH":
        command += ["-n", "IGH_VDJ"]

    PATH_TO_DATA_DIR.mkdir(parents=True, exist_ok=True)

    # the tool picks its own file name, so download into a temporary directory and rename
    with tempfile.TemporaryDirectory() as tmp_dir:
        subprocess.run(command, cwd=tmp_dir, check=True)
        downloaded = list(Path(tmp_dir).glob("*.json"))
        assert len(downloaded) == 1, f"Expected one downloaded file for {locus}, found {len(downloaded)}"
        shutil.move(downloaded[0], path)


def get_ig_aa_sequence_data(paths: Tuple[Path]) -> dict:
    alleles = get_alleles(paths)

    v_gene_sequence_data = get_v_gene_sequence_data(alleles)
    d_gene_sequence_data = get_d_gene_sequence_data(alleles)
    j_gene_sequence_data = get_j_gene_sequence_data(alleles)
    return {**v_gene_sequence_data, **d_gene_sequence_data, **j_gene_sequence_data}


def get_alleles(paths: Tuple[Path]) -> dict:
    alleles = dict()

    for path in paths:
        with open(path, "r") as f:
            germline_sets = json.load(f)["GermlineSet"]

        assert len(germline_sets) == 1, f"Found multiple germline sets in {path}"

        for allele in germline_sets[0]["allele_descriptions"]:
            alleles[allele["label"]] = allele

    return alleles


def get_alleles_of_type(alleles: dict, sequence_type: str) -> dict:
    alleles_of_type = dict()

    for name, allele in alleles.items():
        if allele["sequence_type"] != sequence_type:
            continue

        if not allele["coding_sequence"]:
            print(f"Warning: {name} has no coding sequence, skipping")
            continue

        alleles_of_type[name] = allele

    return alleles_of_type


def get_v_gene_sequence_data(alleles: dict) -> dict:
    aa_seqs = collections.defaultdict(dict)

    for name, allele in get_alleles_of_type(alleles, "V").items():
        nt_seq = allele["coding_sequence"].upper()
        v_region = translate(nt_seq)

        aa_seqs[name]["functionality"] = get_functionality(allele, v_region)
        aa_seqs[name].update(get_v_regions(name, allele, nt_seq))
        aa_seqs[name]["V-REGION"] = v_region

    return script_utility.omit_incomplete_regions(aa_seqs)


def get_d_gene_sequence_data(alleles: dict) -> dict:
    aa_seqs = collections.defaultdict(dict)

    for name, allele in get_alleles_of_type(alleles, "D").items():
        # the D reading frame is arbitrary, so stop codons do not indicate a pseudogene
        aa_seqs[name]["functionality"] = get_functionality(allele)
        aa_seqs[name]["D-REGION"] = translate(allele["coding_sequence"].upper())

    return aa_seqs


def get_j_gene_sequence_data(alleles: dict) -> dict:
    aa_seqs = collections.defaultdict(dict)

    for name, allele in get_alleles_of_type(alleles, "J").items():
        frame = allele["j_codon_frame"]

        if frame is None:
            print(f"Warning: {name} has no J codon frame, skipping")
            continue

        j_region = translate(allele["coding_sequence"].upper()[frame - 1:])

        aa_seqs[name]["functionality"] = get_functionality(allele, j_region)
        aa_seqs[name]["J-REGION"] = j_region

        conserved_aa = get_conserved_aa(name, allele, j_region)

        if conserved_aa == "F":
            aa_seqs[name]["J-PHE"] = conserved_aa
        elif conserved_aa == "W":
            aa_seqs[name]["J-TRP"] = conserved_aa
        elif conserved_aa is not None:
            aa_seqs[name]["J-conserved"] = conserved_aa

    return script_utility.add_j_motifs(aa_seqs)


def get_v_regions(name: str, allele: dict, nt_seq: str) -> dict:
    delineations = [d for d in allele["v_gene_delineations"] or [] if d["delineation_scheme"] == "IMGT"]

    if len(delineations) == 0:
        print(f"Warning: {name} has no IMGT delineation, only V-REGION available")
        return dict()

    delineation = delineations[0]
    regions = dict()

    for stem, label in V_REGIONS:
        start, end = delineation[f"{stem}_start"], delineation[f"{stem}_end"]

        if start is None or end is None:
            print(f"Warning: {name} has no {label} coordinates, region omitted")
            continue

        if (start - 1) % 3 != 0 or (end - start + 1) % 3 != 0:
            print(f"Warning: {name} {label} is not codon aligned, region omitted")
            continue

        regions[label] = translate(nt_seq[start - 1:end])

    return regions


def get_conserved_aa(name: str, allele: dict, j_region: str) -> str:
    # j_cdr3_end marks the first nucleotide of the conserved anchor codon
    if allele["j_cdr3_end"] is None:
        print(f"Warning: {name} has no j_cdr3_end, conserved residue unknown")
        return None

    offset = allele["j_cdr3_end"] - allele["j_codon_frame"]

    if offset < 0 or offset % 3 != 0 or offset // 3 >= len(j_region):
        print(f"Warning: {name} j_cdr3_end is not codon aligned, conserved residue unknown")
        return None

    return j_region[offset // 3]


def get_functionality(allele: dict, aa_seq: str = None) -> str:
    # OGRDB only records functional or not; a stop codon separates P from ORF
    if allele["functional"]:
        return "F"

    if aa_seq is not None and "*" in aa_seq:
        return "P"

    return "ORF"


def translate(nt_seq: str) -> str:
    codons = [nt_seq[i:i + 3] for i in range(0, len(nt_seq) - len(nt_seq) % 3, 3)]
    return "".join(CODON_TABLE[codon] if codon in CODON_TABLE else "X" for codon in codons)


if __name__ == "__main__":
    main()