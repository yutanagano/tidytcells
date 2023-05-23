import json
from importlib.resources import files


with files(__name__).joinpath("homosapiens_tcr.json").open("rb") as f:
    HOMOSAPIENS_TCR = json.load(f)
with files(__name__).joinpath("homosapiens_tcr_synonyms.json").open("rb") as f:
    HOMOSAPIENS_TCR_SYNONYMS = json.load(f)
with files(__name__).joinpath("homosapiens_tcr_aa_sequences.json").open("rb") as f:
    HOMOSAPIENS_TCR_AA_SEQUENCES = json.load(f)
with files(__name__).joinpath("homosapiens_mhc.json").open("rb") as f:
    HOMOSAPIENS_MHC = json.load(f)
with files(__name__).joinpath("homosapiens_mhc_synonyms.json").open("rb") as f:
    HOMOSAPIENS_MHC_SYNONYMS = json.load(f)


with files(__name__).joinpath("musmusculus_tcr.json").open("rb") as f:
    MUSMUSCULUS_TCR = json.load(f)
with files(__name__).joinpath("musmusculus_mhc.json").open("rb") as f:
    MUSMUSCULUS_MHC = json.load(f)
with files(__name__).joinpath("musmusculus_mhc_synonyms.json").open("rb") as f:
    MUSMUSCULUS_MHC_SYNONYMS = json.load(f)


AMINO_ACIDS = frozenset(
    (
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
    )
)
