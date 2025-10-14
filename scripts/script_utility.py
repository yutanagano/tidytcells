from bs4 import BeautifulSoup
import collections
import itertools
from io import StringIO
import json
import pandas as pd
from pandas import DataFrame
from pathlib import Path
import requests
from typing import Tuple
import re
import time

def get_tr_alleles_list(species: str) -> dict:
    gene_groups = (
        "TRAC",
        "TRAV",
        "TRAJ",
        "TRBC",
        "TRBV",
        "TRBD",
        "TRBJ",
        "TRGC",
        "TRGV",
        "TRGJ",
        "TRDC",
        "TRDV",
        "TRDD",
        "TRDJ",
    )

    return get_alleles_list(species, gene_groups)


def get_ig_alleles_list(species: str) -> dict:
    gene_groups = (
        "IGHC",
        "IGHV",
        "IGHD",
        "IGHJ",
        "IGLC",
        "IGLV",
        "IGLJ",
        "IGKC",
        "IGKV",
        "IGKJ",
    )

    return get_alleles_list(species, gene_groups)


def get_alleles_list(species: str, gene_groups: Tuple[str]) -> dict:
    alleles_per_gene_group = [
        get_alleles_for_gene_group_for_species(gene_group, species)
        for gene_group in gene_groups
    ]

    combined = dict()
    for alleles_dict in alleles_per_gene_group:
        combined.update(alleles_dict)

    return combined


def get_alleles_for_gene_group_for_species(gene_group: str, species: str) -> dict:
    alleles = collections.defaultdict(dict)

    time.sleep(5) # to prevent getting blocked
    response = requests.get(
        f"https://www.imgt.org/genedb/GENElect?query=7.2+{gene_group}&species={species}"
    )
    parser = BeautifulSoup(response.text, features="html.parser")
    fasta = parser.find_all("pre")[1].string
    header_lines = filter(lambda line: line.startswith(">"), fasta.splitlines())

    for line in header_lines:
        gene, allele_designation, functionality = parse_fasta_header(line)
        alleles[gene][allele_designation] = functionality

    return alleles


def get_tr_aa_sequence_data(species: str) -> dict:
    v_gene_sequence_data = get_tr_v_gene_sequence_data(species)
    d_gene_sequence_data = get_tr_d_gene_sequence_data(species)
    j_gene_sequence_data = get_tr_j_gene_sequence_data(species)
    return {**v_gene_sequence_data, **d_gene_sequence_data, **j_gene_sequence_data}


def get_ig_aa_sequence_data(species: str) -> dict:
    v_gene_sequence_data = get_ig_v_gene_sequence_data(species)
    d_gene_sequence_data = get_ig_d_gene_sequence_data(species)
    j_gene_sequence_data = get_ig_j_gene_sequence_data(species)
    return {**v_gene_sequence_data, **d_gene_sequence_data, **j_gene_sequence_data}


def get_tr_v_gene_sequence_data(species: str) -> dict:
    gene_groups = ("TRAV", "TRBV", "TRGV", "TRDV")
    return get_v_gene_sequence_data(species, gene_groups)


def get_tr_d_gene_sequence_data(species: str) -> dict:
    gene_groups = ("TRBD", "TRDD")
    return get_d_gene_sequence_data(species, gene_groups)


def get_tr_j_gene_sequence_data(species: str) -> dict:
    gene_groups = ("TRAJ", "TRBJ", "TRGJ", "TRDJ")
    return get_j_gene_sequence_data(species, gene_groups)


def get_ig_v_gene_sequence_data(species: str) -> dict:
    gene_groups = ("IGHV", "IGLV", "IGKV")
    return get_v_gene_sequence_data(species, gene_groups)


def get_ig_d_gene_sequence_data(species: str) -> dict:
    gene_groups = ("IGHD",)
    return get_d_gene_sequence_data(species, gene_groups)


def get_ig_j_gene_sequence_data(species: str) -> dict:
    gene_groups = ("IGHJ", "IGLJ", "IGKJ")
    return get_j_gene_sequence_data(species, gene_groups)


def get_v_gene_sequence_data(species: str, gene_groups: Tuple[str]) -> dict:
    labels = ("FR1-IMGT", "CDR1-IMGT", "FR2-IMGT", "CDR2-IMGT", "FR3-IMGT", "V-REGION")
    return add_v_motifs(get_gene_sequence_data(labels, gene_groups, species))


def get_d_gene_sequence_data(species: str, gene_groups: Tuple[str]) -> dict:
    labels = ("D-REGION",)
    return get_gene_sequence_data(labels, gene_groups, species)


def get_j_gene_sequence_data(species: str, gene_groups: Tuple[str]) -> dict:
    labels = ("J-REGION", "J-MOTIF", "J-PHE", "J-TRP")
    data = get_gene_sequence_data(labels, gene_groups, species)

    # Add known non-canonical cases
    if species == "Homo+sapiens":
        if "TRA" in gene_groups:
            data["TRAJ35*01"]["J-CYS"] = "C"
        if "TRB" in gene_groups:
            data["TRBJ2-7*02"]["J-VAL"] = "V"

    return add_j_motifs(data, species)


def get_gene_sequence_data(
    labels: Tuple[str], gene_groups: Tuple[str], species: str
) -> dict:
    data_per_gene_group_per_label = [
        get_sequence_data_for_label_for_gene_group_for_species(
            label, gene_group, species
        )
        for label, gene_group in itertools.product(labels, gene_groups)
    ]

    combined = collections.defaultdict(dict)
    for alleles_dict in data_per_gene_group_per_label:
        for allele, data in alleles_dict.items():
            combined[allele].update(data)

    return combined


def get_sequence_data_for_label_for_gene_group_for_species(
    label: str, gene_group: str, species: str
) -> dict:
    aa_seqs = collections.defaultdict(dict)

    time.sleep(5)  # to prevent getting blocked
    response = requests.get(
        f"https://www.imgt.org/genedb/GENElect?query=8.2+{gene_group}&species={species}&IMGTlabel={label}",
        headers={"User-Agent": "Mozilla/5.0"},
    )

    parser = BeautifulSoup(response.text, features="html.parser")
    fasta = parser.find_all("pre")[1].string

    current_allele = None
    for line in fasta.splitlines():
        if line.startswith(">"):
            fields = line.split("|")
            allele = fields[1]
            functionality = fields[3]

            current_allele = allele

            aa_seqs[current_allele]["functionality"] = functionality

            continue

        if current_allele is None:
            continue

        if label == "J-PHE" and line.strip() != "F":
            continue

        if label == "J-TRP" and line.strip() != "W":
            continue

        if not label in aa_seqs[current_allele]:
            aa_seqs[current_allele][label] = line.strip()
        else:
            aa_seqs[current_allele][label] += line.strip()

    return aa_seqs


def parse_fasta_header(line: str) -> Tuple[str]:
    fields = line.split("|")
    allele_name = fields[1]
    gene, allele_designation = allele_name.split("*")
    functionality = fields[3].strip("()[]")

    return gene, allele_designation, functionality


def add_v_motifs(v_aa_dict):
    for allele, seq_data in v_aa_dict.items():
        if "FR3-IMGT" in seq_data:
            if seq_data["FR3-IMGT"].endswith("C") and len(seq_data["FR3-IMGT"]) >= 4:
                v_aa_dict[allele]["V-MOTIF"] = seq_data["FR3-IMGT"][-4:]

                if seq_data["FR3-IMGT"] in seq_data["V-REGION"]:
                    cdr3_start_motif = seq_data["V-REGION"][seq_data["V-REGION"].index(seq_data["FR3-IMGT"]) + len(seq_data["FR3-IMGT"]) - 1:]

                    if len(cdr3_start_motif) > 0:
                        v_aa_dict[allele]["V-CDR3-START"] = cdr3_start_motif.rstrip("*")

    return v_aa_dict


def get_motif_idx(j_region, conserved_aa):
    if conserved_aa is not None:
        if conserved_aa not in j_region:
            return None

        if j_region.count(conserved_aa) == 1:
            return j_region.index(conserved_aa)

        if j_region.count(conserved_aa + "G") == 1:  # G is a very common second amino acid in the motif
            return j_region.index(conserved_aa + "G")

    cons_aas_regex = "[FW]" if conserved_aa is None else conserved_aa
    motif_regex = cons_aas_regex + "[AGS][A-Z]G"

    match = re.search(motif_regex, j_region)
    if match:
        return match.start()


def add_j_motifs(j_aa_dict, species):
    for allele, seq_data in j_aa_dict.items():
        if "J-REGION" not in seq_data or len(seq_data["J-REGION"]) < 4:
            continue

        conserved_aa = seq_data["J-PHE"] if "J-PHE" in seq_data \
            else seq_data["J-TRP"] if "J-TRP" in seq_data \
            else seq_data["J-CYS"] if "J-CYS" in seq_data \
            else seq_data["J-VAL"] if "J-VAL" in seq_data \
            else None

        if not ("J-MOTIF" in seq_data and len(seq_data["J-MOTIF"]) == 4):
            motif_idx = get_motif_idx(seq_data["J-REGION"], conserved_aa)

            if motif_idx is not None:
                if "J-MOTIF" in seq_data:
                    if seq_data["J-REGION"][motif_idx: motif_idx + 4] != seq_data["J-MOTIF"]:
                        pass
                seq_data["J-MOTIF"] = seq_data["J-REGION"][motif_idx: motif_idx + 4]

        if "J-MOTIF" in seq_data and len(seq_data["J-MOTIF"]) == 4:
            if conserved_aa is None:
                conserved_aa = seq_data["J-MOTIF"][0]
                if conserved_aa == "F":
                    seq_data["J-PHE"] = "F"
                elif conserved_aa == "W":
                    seq_data["J-TRP"] = "W"

            cdr3_end_motif = seq_data["J-REGION"][0:seq_data["J-REGION"].index(seq_data["J-MOTIF"]) + 1]

            if len(cdr3_end_motif) > 1:
                j_aa_dict[allele]["J-CDR3-END"] = cdr3_end_motif.lstrip("*")

    return j_aa_dict


def save_as_json(data: dict, file_name: str) -> None:
    path_to_resources_dir = Path("src") / "tidytcells" / "_resources"
    with open(path_to_resources_dir / file_name, "w") as f:
        json.dump(data, f, indent=4)


def fetch_hgnc_data() -> DataFrame:
    response = requests.get(
        "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_locus_type&col=family.name&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit",
        stream=True,
    )
    buffer = StringIO(response.text)
    return pd.read_csv(buffer, sep="\t", on_bad_lines="warn")
