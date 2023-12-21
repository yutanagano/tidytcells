from bs4 import BeautifulSoup
import collections
import json
from pathlib import Path
import requests


def get_tr_alleles_list(species: str) -> dict:
    gene_groups = (
        "TRAC", "TRAJ", "TRAV",
        "TRBC", "TRBJ", "TRBD", "TRBV",
        "TRGC", "TRGJ", "TRGV",
        "TRDC", "TRDJ", "TRDD", "TRDV"
    )

    alleles_per_gene_group = [get_tr_alleles_for_gene_group_for_species(gene_group, species) for gene_group in gene_groups]

    combined = dict()
    for alleles_dict in alleles_per_gene_group:
        combined.update(alleles_dict)
    
    return combined


def get_tr_alleles_for_gene_group_for_species(gene_group: str, species: str) -> dict:
    alleles = collections.defaultdict(dict)

    response = requests.get(f"https://www.imgt.org/genedb/GENElect?query=7.2+{gene_group}&species={species}")
    parser = BeautifulSoup(response.text, features="html.parser")
    fasta = parser.find_all("pre")[1].string
    header_lines = filter(lambda line: line.startswith(">"), fasta.splitlines())

    for line in header_lines:
        gene, allele_designation, functionality = parse_fasta_header(line)
        alleles[gene][allele_designation] = functionality
    
    return alleles


def parse_fasta_header(line: str) -> (str, str, str):
    fields = line.split("|")
    allele_name = fields[1]
    gene, allele_designation = allele_name.split("*")
    functionality = fields[3].strip("()[]")

    return gene, allele_designation, functionality


def save_as_json(data: dict, file_name: str) -> None:
    path_to_resources_dir = Path("src")/"tidytcells"/"_resources"
    with open(path_to_resources_dir/file_name, "w") as f:
        json.dump(data, f, indent=4)