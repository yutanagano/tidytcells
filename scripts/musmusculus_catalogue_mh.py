from pathlib import Path
import json
import pandas as pd
import re


df = pd.read_excel(Path("data") / "musmusculus_mh.ods").drop(
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


# ## List all MH genes


# Get all genes
genes = df["gene name"].unique().tolist()


with open(
    Path("src") / "tidytcells" / "_resources" / "valid_musmusculus_mh.json", "w"
) as f:
    json.dump(genes, f, indent=4)


# ## Get deprecated names/synonyms


mh_synonyms = df[["gene name", "synonym"]]

# Group together by synonym
mh_synonyms = mh_synonyms.groupby("synonym").aggregate(lambda x: x.tolist())


# Discard ambiguous synonyms
mh_synonyms = mh_synonyms[mh_synonyms["gene name"].map(len) == 1].copy()
mh_synonyms["gene name"] = mh_synonyms["gene name"].map(lambda x: x.pop())


# Discard redundant items (synonym == approved symbol)
mh_synonyms = mh_synonyms[mh_synonyms.index != mh_synonyms["gene name"]]


# Discard synonyms that are also names of other valid genes
mh_synonyms = mh_synonyms[mh_synonyms.index.map(lambda x: x not in genes)]


mh_synonyms["gene name"].to_json(
    Path("src") / "tidytcells" / "_resources" / "musmusculus_mh_synonyms.json",
    indent=4,
)
