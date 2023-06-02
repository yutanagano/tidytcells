from importlib import resources
import json
import sys

if sys.version_info >= (3, 9):

    def get_json_resource(filename: str) -> dict:
        with resources.files(__name__).joinpath(filename).open("r") as f:
            return json.load(f)

else:

    def get_json_resource(filename: str) -> dict:
        with resources.open_text(__name__, filename) as f:
            return json.load(f)


HOMOSAPIENS_TCR = get_json_resource("homosapiens_tcr.json")
HOMOSAPIENS_TCR_SYNONYMS = get_json_resource("homosapiens_tcr_synonyms.json")
HOMOSAPIENS_TCR_AA_SEQUENCES = get_json_resource("homosapiens_tcr_aa_sequences.json")
HOMOSAPIENS_MHC = get_json_resource("homosapiens_mhc.json")
HOMOSAPIENS_MHC_SYNONYMS = get_json_resource("homosapiens_mhc_synonyms.json")


MUSMUSCULUS_TCR = get_json_resource("musmusculus_tcr.json")
MUSMUSCULUS_MHC = get_json_resource("musmusculus_mhc.json")
MUSMUSCULUS_MHC_SYNONYMS = get_json_resource("musmusculus_mhc_synonyms.json")


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
