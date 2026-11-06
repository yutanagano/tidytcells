from pathlib import Path
from io import StringIO
import pandas as pd
import requests
import re
from numpy import ndarray
import script_utility
from typing import Iterable


def main() -> None:
    print("Fetching HLA allele data from IMGT...")
    valid_allele_tree = get_hla_alleles_tree()
    script_utility.save_as_json(valid_allele_tree, "ipd_homosapiens_mh_valid_allele.json")

    print("Fetching HLA symbol synonyms from HGNC...")
    synonyms_data = script_utility.get_hgnc_synonyms_data(valid_allele_tree)
    script_utility.save_as_json(synonyms_data, "ipd_homosapiens_mh_synonyms_hgnc.json")


def get_hla_alleles_tree() -> dict:
    alleles = get_alleles_data()
    g_groups = get_g_groups_data()
    p_groups = get_p_groups_data()

    alleles_decomposed = [decompose_hla(allele) for allele in alleles]
    g_groups_decomposed = [decompose_hla(group, 4) for group in g_groups]
    p_groups_decomposed = [decompose_hla(group, 4) for group in p_groups]

    combined_decomposed = list(
        dict.fromkeys(
            sorted(alleles_decomposed + g_groups_decomposed + p_groups_decomposed)
        )
    )

    hla_tree = dict()
    make_hla_tree(hla_tree, combined_decomposed)

    return hla_tree


def get_alleles_data() -> ndarray:
    response = requests.get(
        "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom.txt",
        stream=True,
    )
    buffer = StringIO(response.text)
    df = pd.read_csv(
        buffer,
        sep=";",
        skiprows=6,
        names=[
            "locus",
            "allele",
            "assigned",
            "deleted",
            "identical to",
            "reason for deletion",
        ],
    )
    df = df[df.locus.str.endswith("*")]
    df = df[df.deleted.isna()]
    return df.apply(lambda row: "HLA-" + row["locus"] + row["allele"], axis=1).unique()


def get_g_groups_data() -> ndarray:
    response = requests.get(
        "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom_g.txt",
        stream=True,
    )
    buffer = StringIO(response.text)
    df = pd.read_csv(buffer, sep=";", skiprows=6, names=["locus", "allele", "group"])
    df = df[df.locus.str.endswith("*")]
    df = df[df.group.notna()]
    return df.apply(lambda row: "HLA-" + row["locus"] + row["group"], axis=1).unique()


def get_p_groups_data() -> ndarray:
    response = requests.get(
        "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom_p.txt",
        stream=True,
    )
    buffer = StringIO(response.text)
    df = pd.read_csv(buffer, sep=";", skiprows=6, names=["locus", "allele", "group"])
    df = df[df.locus.str.endswith("*")]
    df = df[df.group.notna()]
    return df.apply(lambda row: "HLA-" + row["locus"] + row["group"], axis=1).unique()


def decompose_hla(gene_str: str, max_spec_field_depth: int = 2):
    m = re.match(r"^([A-Za-z0-9\-]+)\*([\dGP:]+)[LSCAQN]?$", gene_str)

    if m is None:
        raise ValueError(gene_str)

    gene = m.group(1)
    spec_fields = m.group(2).split(":")[:max_spec_field_depth]

    return (gene,) + tuple(spec_fields)


def make_hla_tree(current_root: dict, token_lists: list) -> None:
    first_tokens = list(
        dict.fromkeys(sorted([token_list[0] for token_list in token_lists]))
    )

    for token in first_tokens:
        current_root[token] = {}

        new_token_lists = [
            token_list[1:]
            for token_list in token_lists
            if token_list[0] == token and len(token_list) > 1
        ]

        make_hla_tree(current_root[token], new_token_lists)



if __name__ == "__main__":
    main()
