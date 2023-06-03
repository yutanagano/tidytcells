from pathlib import Path
import json
import pandas as pd


alleles = dict()

with open(Path("data") / "homosapiens_tcr.fasta", "r") as f:
    for line in f.readlines():
        if line.startswith(">"):
            fields = line.split("|")
            allele_name = fields[1]
            gene = allele_name.split("*")[0]
            allele_designation = allele_name.split("*")[1]
            functionality = fields[3].strip("()[]")

            if not gene in alleles:
                alleles[gene] = dict()

            alleles[gene][allele_designation] = functionality


with open(Path("src") / "tidytcells" / "_resources" / "homosapiens_tcr.json", "w") as f:
    json.dump(alleles, f, indent=4)


# ## Get deprecated names/synonyms


hgnc = pd.read_csv(Path("data") / "hgnc.tsv", sep="\t")


# Get only TCR genes
tcr_genes = hgnc[hgnc["Locus type"].str.contains("T cell receptor gene")].copy()


# Put back slashes behind DV designations and OR designations
tcr_genes["Approved symbol"] = tcr_genes["Approved symbol"].str.replace(
    r"(?<!TR)DV", "/DV", regex=True
)
tcr_genes["Approved symbol"] = tcr_genes["Approved symbol"].str.replace(
    r"OR", "/OR", regex=True
)


# Only keep genes whose 'approved symbols' are in our IMGT list
tcr_genes = tcr_genes[tcr_genes["Approved symbol"].map(lambda x: x in alleles)].copy()


# Get TCR genes with aliases
tcr_genes_with_aliases = tcr_genes[tcr_genes["Alias symbols"].notna()][
    ["Approved symbol", "Alias symbols"]
]
tcr_genes_with_aliases["Alias symbols"] = tcr_genes_with_aliases["Alias symbols"].map(
    lambda x: x.split(", ")
)
tcr_genes_with_aliases.columns = ["Approved symbol", "Synonym"]
tcr_genes_with_aliases = tcr_genes_with_aliases.explode("Synonym")


# Get TCR genes with deprecated names
tcr_genes_with_depnames = tcr_genes[tcr_genes["Previous symbols"].notna()][
    ["Approved symbol", "Previous symbols"]
]
tcr_genes_with_depnames["Previous symbols"] = tcr_genes_with_depnames[
    "Previous symbols"
].map(lambda x: x.split(", "))
tcr_genes_with_depnames.columns = ["Approved symbol", "Synonym"]
tcr_genes_with_depnames = tcr_genes_with_depnames.explode("Synonym")


# Combine both tables
tcr_synonyms = pd.concat([tcr_genes_with_aliases, tcr_genes_with_depnames])

# Remove any names that are now redundant (approved symbol and synonym are the same)
tcr_synonyms = tcr_synonyms[tcr_synonyms["Approved symbol"] != tcr_synonyms["Synonym"]]

# Group together by synonym
tcr_synonyms = tcr_synonyms.groupby("Synonym").aggregate(lambda x: x.tolist())


# Remove ambiguous synonyms
tcr_synonyms = tcr_synonyms[tcr_synonyms["Approved symbol"].map(len) == 1].copy()
tcr_synonyms["Approved symbol"] = tcr_synonyms["Approved symbol"].map(lambda x: x.pop())
tcr_synonyms.index = tcr_synonyms.index.str.upper()


# Remove any synonyms that are also names of other valid genes
tcr_synonyms = tcr_synonyms[tcr_synonyms.index.map(lambda x: x not in alleles)]


tcr_synonyms["Approved symbol"].to_json(
    Path("src") / "tidytcells" / "_resources" / "homosapiens_tcr_synonyms.json",
    indent=4,
)


# ## Catalogue amino acid sequences for TCR V genes


v_aa_seqs = dict()

with open(Path("data") / "homosapiens_tcr_aa_sequences.fasta", "r") as f:
    current_v = None
    current_segment = None

    for line in f.readlines():
        if line.startswith(">"):
            fields = line.split("|")
            current_v = fields[1]
            current_segment = fields[4]
            continue

        if not current_v in v_aa_seqs:
            v_aa_seqs[current_v] = dict()

        if not current_segment in v_aa_seqs[current_v]:
            v_aa_seqs[current_v][current_segment] = line.strip()
        else:
            v_aa_seqs[current_v][current_segment] += line.strip()


with open(
    Path("src") / "tidytcells" / "_resources" / "homosapiens_tcr_aa_sequences.json", "w"
) as f:
    json.dump(v_aa_seqs, f, indent=4)
