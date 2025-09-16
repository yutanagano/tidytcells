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

    if species == "Homo+sapiens" and "TRA" in gene_groups:
        data["TRAJ35*01"]["J-CYS"] = "C"

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

            if "F" in functionality:
                current_allele = allele
            else:
                current_allele = None

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
                    _v, cdr3_motif = seq_data["V-REGION"].split(seq_data["FR3-IMGT"])
                    if len(cdr3_motif) > 0:
                        v_aa_dict[allele]["V-CDR3-MOTIF"] = cdr3_motif

    return v_aa_dict


def add_j_motifs(j_aa_dict, species):
    for allele, seq_data in j_aa_dict.items():
        if "J-REGION" not in seq_data or len(seq_data["J-REGION"]) < 4:
            continue

        if not ("J-MOTIF" in seq_data and len(seq_data["J-MOTIF"]) == 4):
            conserved_aa = seq_data["J-PHE"] if "J-PHE" in seq_data \
                else seq_data["J-TRP"] if "J-TRP" in seq_data \
                else seq_data["J-CYS"] if "J-CYS" in seq_data \
                else None

            if conserved_aa is None or conserved_aa not in seq_data["J-REGION"]:
                continue

            if seq_data["J-REGION"].count(conserved_aa) == 1:
                motif_idx = seq_data["J-REGION"].index(conserved_aa)
            elif seq_data["J-REGION"].count(conserved_aa + "G") == 1: # G is a very common second amino acid in the motif
                motif_idx = seq_data["J-REGION"].index(conserved_aa + "G")
            elif seq_data["J-REGION"][1:].count(conserved_aa) == 1:  # J-REGION sometimes starts with F, which is unlikely to be the motif-F
                motif_idx = seq_data["J-REGION"][1:].index(conserved_aa) + 1
            elif species == "Mus+musculus" and allele.startswith("TRG") and seq_data["J-REGION"].count(conserved_aa + "A") == 1:
                motif_idx = seq_data["J-REGION"].index(conserved_aa + "A")
            else:
                continue

            seq_data["J-MOTIF"] = seq_data["J-REGION"][motif_idx: motif_idx + 4]

        cdr3_motif, _j = seq_data["J-REGION"].split(seq_data["J-MOTIF"])

        if len(cdr3_motif) > 0:
            j_aa_dict[allele]["J-CDR3-MOTIF"] = cdr3_motif

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
