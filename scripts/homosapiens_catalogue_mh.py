from pathlib import Path
import json
import pandas as pd
import re
from urllib.request import urlopen


PROJECT_PATH = Path.cwd().parent.resolve()

with urlopen(
    "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom.txt"
) as f:
    alleles_df = pd.read_csv(
        f,
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
    alleles_df = alleles_df[alleles_df["locus"].str.endswith("*")]

with urlopen(
    "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom_g.txt"
) as f:
    g_df = pd.read_csv(f, sep=";", skiprows=6, names=["locus", "allele", "group"])
    g_df = g_df[g_df["locus"].str.endswith("*")]

with urlopen(
    "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom_p.txt"
) as f:
    p_df = pd.read_csv(f, sep=";", skiprows=6, names=["locus", "allele", "group"])
    p_df = p_df[p_df["locus"].str.endswith("*")]

alleles = alleles_df[alleles_df["deleted"].isna()]
alleles = alleles.apply(
    lambda row: "HLA-" + row["locus"] + row["allele"], axis=1
).unique()

g_groups = g_df[g_df["group"].notna()]
g_groups = g_groups.apply(
    lambda row: "HLA-" + row["locus"] + row["group"], axis=1
).unique()

p_groups = p_df[p_df["group"].notna()]
p_groups = p_groups.apply(
    lambda row: "HLA-" + row["locus"] + row["group"], axis=1
).unique()


def decompose_hla(gene_str: str, max_spec_field_depth: int = 2):
    m = re.match(r"^([A-Z0-9\-]+)\*([\dGP:]+)[LSCAQN]?$", gene_str)

    if m is None:
        raise ValueError(gene_str)

    gene = m.group(1)
    spec_fields = m.group(2).split(":")[:max_spec_field_depth]

    return (gene,) + tuple(spec_fields)


g_groups_decomposed = [
    decompose_hla(g_group, 4)
    for g_group in g_groups
    if decompose_hla(g_group, 4) is not None
]
p_groups_decomposed = [
    decompose_hla(p_group, 4)
    for p_group in p_groups
    if decompose_hla(p_group, 4) is not None
]
proteins_decomposed = [decompose_hla(allele) for allele in alleles]

combined_decomposed = list(
    dict.fromkeys(
        sorted(g_groups_decomposed + p_groups_decomposed + proteins_decomposed)
    )
)


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


hla_tree = {}

make_hla_tree(hla_tree, combined_decomposed)

with open(
    Path("src") / "tidytcells" / "_resources" / "valid_homosapiens_mh.json", "w"
) as f:
    json.dump(hla_tree, f, indent=4)


hgnc = pd.read_csv(Path("data") / "hgnc.tsv", sep="\t")

# Get only MH genes
mh_genes = hgnc[hgnc["Gene group name"].notna()]
mh_genes = mh_genes[
    mh_genes["Gene group name"].str.contains("Histocompatibility complex")
]

# Only keep genes whose 'approved symbols' are in our IMGT list
mh_genes = mh_genes[mh_genes["Approved symbol"].map(lambda x: x in hla_tree)]

# Get MH genes with aliases
mh_genes_with_aliases = mh_genes[mh_genes["Alias symbols"].notna()][
    ["Approved symbol", "Alias symbols"]
]
mh_genes_with_aliases["Alias symbols"] = mh_genes_with_aliases["Alias symbols"].map(
    lambda x: x.split(", ")
)
mh_genes_with_aliases.columns = ["Approved symbol", "Synonym"]
mh_genes_with_aliases = mh_genes_with_aliases.explode("Synonym")

# Get MH genes with deprecated names
mh_genes_with_depnames = mh_genes[mh_genes["Previous symbols"].notna()][
    ["Approved symbol", "Previous symbols"]
]
mh_genes_with_depnames["Previous symbols"] = mh_genes_with_depnames[
    "Previous symbols"
].map(lambda x: x.split(", "))
mh_genes_with_depnames.columns = ["Approved symbol", "Synonym"]
mh_genes_with_depnames = mh_genes_with_depnames.explode("Synonym")

# Combine both tables
mh_synonyms = pd.concat([mh_genes_with_aliases, mh_genes_with_depnames])

# Capitalise synonyms
mh_synonyms["Synonym"] = mh_synonyms["Synonym"].str.upper()

# Group together by synonym
mh_synonyms = mh_synonyms.groupby("Synonym").aggregate(lambda x: x.tolist())

# Discard ambiguous synonyms
mh_synonyms = mh_synonyms[mh_synonyms["Approved symbol"].map(len) == 1].copy()
mh_synonyms["Approved symbol"] = mh_synonyms["Approved symbol"].map(lambda x: x.pop())

# Discard redundant items (synonym == approved symbol)
mh_synonyms = mh_synonyms[mh_synonyms.index != mh_synonyms["Approved symbol"]]

# Discard synonyms that are also names of other valid genes
mh_synonyms = mh_synonyms[mh_synonyms.index.map(lambda x: x not in hla_tree)]

mh_synonyms["Approved symbol"].to_json(
    Path("src") / "tidytcells" / "_resources" / "homosapiens_mh_synonyms.json",
    indent=4,
)
