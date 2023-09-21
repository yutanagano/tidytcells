from pathlib import Path
import json
import pandas as pd
import re


df = pd.read_excel(Path("data") / "musmusculus_mhc.ods").drop(
    columns=["Unnamed: 2", "Unnamed: 4"]
)
df.columns = ["group", "subgroup", "gene name", "synonym"]
df["group"] = df["group"].ffill()
df["subgroup"] = df["subgroup"].ffill()

df["synonym"] = df["synonym"].map(lambda s: re.sub(r"[\[\(].*[\]\)]", "", s))
df["synonym"] = df["synonym"].map(
    lambda s: ["".join(syn.split()).upper().replace("-", "") for syn in s.split(",")]
)
df = df.explode("synonym", ignore_index=True)

df = df.map(lambda x: x.strip())
df = df.drop_duplicates()


df.head()


# ## List all MHC genes


# Get all genes
genes = df["gene name"].unique().tolist()


with open(
    Path("src") / "tidytcells" / "_resources" / "valid_musmusculus_mh.json", "w"
) as f:
    json.dump(genes, f, indent=4)


# ## Get deprecated names/synonyms


mhc_synonyms = df[["gene name", "synonym"]]

# Group together by synonym
mhc_synonyms = mhc_synonyms.groupby("synonym").aggregate(lambda x: x.tolist())


# Discard ambiguous synonyms
mhc_synonyms = mhc_synonyms[mhc_synonyms["gene name"].map(len) == 1].copy()
mhc_synonyms["gene name"] = mhc_synonyms["gene name"].map(lambda x: x.pop())


# Discard redundant items (synonym == approved symbol)
mhc_synonyms = mhc_synonyms[mhc_synonyms.index != mhc_synonyms["gene name"]]


# Discard synonyms that are also names of other valid genes
mhc_synonyms = mhc_synonyms[mhc_synonyms.index.map(lambda x: x not in genes)]


mhc_synonyms["gene name"].to_json(
    Path("src") / "tidytcells" / "_resources" / "musmusculus_mh_synonyms.json",
    indent=4,
)
