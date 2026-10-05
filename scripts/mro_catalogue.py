import collections
import re
import shutil
import tempfile
import urllib.request
import zipfile
from pathlib import Path

import pandas as pd

import script_utility


PATH_TO_DATA_DIR = Path("data") / "MRO"

RELEASE_URL = "https://github.com/IEDB/MRO/releases/latest/download"

SPECIES_TAXON = {"homosapiens": "NCBITaxon:9606"}

ALLELE_PATTERN = re.compile(r"^([A-Za-z0-9]+-[A-Z0-9]+)\*(\d+(?::\d+)*)[NLSCAQ]?(?: .+ mutant)?$")


def main() -> None:
    print("Fetching MRO data...")
    protein_complexes_path = get_protein_complexes_path()
    molecule_export_path = get_molecule_export_path()

    for species in SPECIES_TAXON:
        print(f"Converting MRO data for {species}...")
        df = load_molecules(protein_complexes_path, molecule_export_path, species)

        script_utility.save_as_json(build_valid_tree(df["IEDB alternative term"]), f"valid_{species}_mh_mro.json")
        script_utility.save_as_json(build_synonyms(df), f"{species}_mh_synonyms_allele_mro.json")


def get_protein_complexes_path() -> Path:
    path = PATH_TO_DATA_DIR / "mro-protein-complexes.tsv"

    if not path.exists():
        PATH_TO_DATA_DIR.mkdir(parents=True, exist_ok=True)
        download(f"{RELEASE_URL}/mro-protein-complexes.tsv", path)

    return path


def get_molecule_export_path() -> Path:
    path = PATH_TO_DATA_DIR / "molecule_export.tsv"

    if not path.exists():
        PATH_TO_DATA_DIR.mkdir(parents=True, exist_ok=True)

        with tempfile.TemporaryDirectory() as tmp_dir:
            zip_path = Path(tmp_dir) / "iedb.zip"
            download(f"{RELEASE_URL}/iedb.zip", zip_path)

            with zipfile.ZipFile(zip_path) as archive:
                members = [m for m in archive.namelist() if Path(m).name == path.name]
                assert len(members) == 1, f"Expected one {path.name} in iedb.zip, found {len(members)}"

                with archive.open(members[0]) as source, open(path, "wb") as target:
                    shutil.copyfileobj(source, target)

    return path


def download(url: str, path: Path) -> None:
    with urllib.request.urlopen(url) as response, open(path, "wb") as f:
        shutil.copyfileobj(response, f)


def load_molecules(protein_complexes_path: Path, molecule_export_path: Path, species: str) -> pd.DataFrame:
    complexes = pd.read_csv(protein_complexes_path, sep="\t", dtype=str).fillna("")
    molecules = pd.read_csv(molecule_export_path, sep="\t", dtype=str).fillna("")
    df = complexes.merge(molecules[["MRO ID", "In Taxon ID"]], left_on="id", right_on="MRO ID")
    return df[df["In Taxon ID"] == SPECIES_TAXON[species]]


def chain_terms(term: str) -> list:
    if "/" not in term:
        return [term]

    first, second = term.split("/", 1)
    prefix = first.split("-", 1)[0]
    return [first, second if second.startswith(f"{prefix}-") else f"{prefix}-{second}"]


def build_valid_tree(terms) -> dict:
    tree = dict()

    for term in terms:
        for chain in chain_terms(term):
            m = ALLELE_PATTERN.match(chain)

            if not m:
                continue

            node = tree.setdefault(m.group(1), dict())
            for field in m.group(2).split(":"):
                node = node.setdefault(field, dict())

    return sort_tree(tree)


def sort_tree(tree: dict) -> dict:
    return {k: sort_tree(tree[k]) for k in sorted(tree)}


def build_synonyms(df: pd.DataFrame) -> dict:
    candidates = collections.defaultdict(set)

    for term, alternatives in zip(df["IEDB alternative term"], df["alternative term"]):
        if "/" in term or "mutant" in term or not ALLELE_PATTERN.match(term):
            continue

        for alternative in filter(None, alternatives.split("|")):
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