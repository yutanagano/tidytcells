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
    script_utility.save_as_json(valid_allele_tree, "valid_homosapiens_mh.json")

    print("Fetching HLA symbol synonyms from HGNC...")
    synonyms_data = get_synonyms_data(valid_allele_tree)
    script_utility.save_as_json(synonyms_data, "homosapiens_mh_synonyms.json")


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


def get_synonyms_data(valid_alleles: Iterable[str]) -> dict:
    hgnc = script_utility.fetch_hgnc_data()

    mh_genes = hgnc[hgnc["Gene group name"].notna()]
    mh_genes = mh_genes[
        mh_genes["Gene group name"].str.contains("Histocompatibility complex")
    ]

    # Only keep genes whose 'approved symbols' are in our IMGT list
    mh_genes = mh_genes[mh_genes["Approved symbol"].map(lambda x: x in valid_alleles)]

    # Get MH genes with "alias symbols"
    mh_genes_with_aliases = mh_genes[mh_genes["Alias symbols"].notna()][
        ["Approved symbol", "Alias symbols"]
    ]
    mh_genes_with_aliases["Alias symbols"] = mh_genes_with_aliases["Alias symbols"].map(
        lambda x: [element.strip() for element in x.split(",")]
    )
    mh_genes_with_aliases.columns = ["Approved symbol", "Synonym"]
    mh_genes_with_aliases = mh_genes_with_aliases.explode("Synonym")

    # Get MH genes with "previous symbols" (deprecated symbols)
    mh_genes_with_depnames = mh_genes[mh_genes["Previous symbols"].notna()][
        ["Approved symbol", "Previous symbols"]
    ]
    mh_genes_with_depnames["Previous symbols"] = mh_genes_with_depnames[
        "Previous symbols"
    ].map(lambda x: [element.strip() for element in x.split(",")])
    mh_genes_with_depnames.columns = ["Approved symbol", "Synonym"]
    mh_genes_with_depnames = mh_genes_with_depnames.explode("Synonym")

    # Combine and uppercase
    mh_synonyms = pd.concat([mh_genes_with_aliases, mh_genes_with_depnames])
    mh_synonyms["Synonym"] = mh_synonyms["Synonym"].str.upper()

    # Remove ambiguous synonyms
    mh_synonyms = mh_synonyms.groupby("Synonym").aggregate(lambda x: x.tolist())
    mh_synonyms = mh_synonyms[mh_synonyms["Approved symbol"].map(len) == 1].copy()
    mh_synonyms["Approved symbol"] = mh_synonyms["Approved symbol"].map(
        lambda x: x.pop()
    )

    # Remove redundant items
    mh_synonyms = mh_synonyms[mh_synonyms.index != mh_synonyms["Approved symbol"]]

    # Remove any synonyms that are also names of other valid genes
    mh_synonyms = mh_synonyms[mh_synonyms.index.map(lambda x: x not in valid_alleles)]

    return mh_synonyms["Approved symbol"].to_dict()


if __name__ == "__main__":
    main()
