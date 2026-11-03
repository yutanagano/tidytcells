import collections
import io
import re
import urllib.request
import zipfile
from pathlib import Path
from typing import List, Optional, Tuple
import pandas as pd

import script_utility


RELEASE_URL = "https://github.com/IEDB/MRO/releases/latest/download"

# Only these species will be included
SPECIES_NAMES = {
    "cattle": "bostaurus",
    "chimpanzee": "pantroglodytes",
    "human": "homosapiens",
    "mouse": "musmusculus",
    "pig": "susscrofa",
    "rat": "rattusnorvegicus",
    "rhesusmacaque": "macacamulatta",
    # "crabeatingmacaque": "macacafascicularis", # to be added, waiting on MRO fix
}


def main() -> None:
    print("Fetching MRO data...")
    molecules = load_molecules(fetch_protein_complexes(), fetch_molecule_export())
    molecules["species"] = [get_species_name(taxon) for taxon in molecules["In Taxon"]]

    for species, df in molecules[molecules["species"].notna()].groupby("species"):
        print(f"Converting MRO data for {species}...")

        df = drop_incomplete_names(df)
        valid = build_valid_symbols(df["IEDB alternative term"])

        if not valid:
            continue

        allele_synonyms = build_mro_synonyms(df)

        script_utility.save_as_json(valid, f"mro_{species}_mh_valid_protein.json")
        script_utility.save_as_json(allele_synonyms, f"mro_{species}_mh_synonyms.json")

        if species == "homosapiens":
            print(f"Adding HGNC gene synonyms for homosapiens...")
            script_utility.save_as_json(script_utility.get_hgnc_synonyms_data(valid), f"mro_{species}_mh_synonyms_hgnc.json")


def fetch_protein_complexes() -> pd.DataFrame:
    return read_tsv(io.BytesIO(download(f"{RELEASE_URL}/mro-protein-complexes.tsv")))


def fetch_molecule_export() -> pd.DataFrame:
    file_name = "molecule_export.tsv"

    with zipfile.ZipFile(io.BytesIO(download(f"{RELEASE_URL}/iedb.zip"))) as archive:
        members = [m for m in archive.namelist() if Path(m).name == file_name]
        assert len(members) == 1, f"Expected one {file_name} in iedb.zip, found {len(members)}"

        with archive.open(members[0]) as f:
            return read_tsv(f)


def download(url: str) -> bytes:
    with urllib.request.urlopen(url) as response:
        return response.read()


def read_tsv(file) -> pd.DataFrame:
    return pd.read_csv(file, sep="\t", dtype=str).fillna("")


def load_molecules(complexes: pd.DataFrame, molecules: pd.DataFrame) -> pd.DataFrame:
    return complexes.merge(molecules[["MRO ID", "In Taxon", "In Taxon ID"]], left_on="id", right_on="MRO ID")


def get_species_name(taxon: str) -> Optional[str]:
    return SPECIES_NAMES.get(re.sub(r"[^a-z]", "", taxon.lower()))


def drop_incomplete_names(df: pd.DataFrame) -> pd.DataFrame:
    # names starting with the organism's common name (e.g. "human CD1a", "mouse") are dropped
    is_symbol = [
        not name.lower().startswith(taxon.lower())
        for name, taxon in zip(df["IEDB alternative term"], df["In Taxon"])
    ]
    return df[is_symbol]


def chain_terms(name: str) -> List[str]:
    if "/" not in name:
        return [name]

    first, second = name.split("/", 1)
    prefix = first.split("-", 1)[0]
    return [first, second if second.startswith(f"{prefix}-") else f"{prefix}-{second}"]


def split_symbol(name: str) -> Tuple[str, List[str]]:
    # drops a mutant description ("HLA-A*02:01 K66A mutant"), and an expression letter ("HLA-A*02:43N"),
    # as tidytcells also removes both from the input before validating it
    symbol = name.split(" ")[0]
    gene, _, designation = symbol.partition("*")
    designation = re.sub(r"(?<=\d)[A-Z]$", "", designation)

    return gene, designation.split(":") if designation else []


def build_valid_symbols(names) -> dict:
    tree = dict()

    for name in names:
        for chain in chain_terms(name):
            gene, fields = split_symbol(chain)

            node = tree.setdefault(gene, dict())
            for field in fields:
                node = node.setdefault(field, dict())

    return sort_tree(tree)


def sort_tree(tree: dict) -> dict:
    return {k: sort_tree(tree[k]) for k in sorted(tree)}


def build_mro_synonyms(df: pd.DataFrame) -> dict:
    candidates = collections.defaultdict(set)

    for term, alternatives in zip(df["IEDB alternative term"], df["alternative term"]):
        if " " in term or "/" in term:
            continue

        for alternative in filter(None, alternatives.split("|")):
            # skips junk alternative names without any digits, e.g. "A" for HLA-DRB3*01:52
            if not re.search(r"\d", alternative):
                continue

            key = alternative.upper()
            if key != term.upper():
                candidates[key].add(term)

    synonyms = dict()

    for key, targets in sorted(candidates.items()):
        if len(targets) > 1:
            print(f"Warning: synonym {key} maps to {sorted(targets)}, skipping")
            continue

        synonyms[key] = targets.pop()

    return synonyms





if __name__ == "__main__":
    main()