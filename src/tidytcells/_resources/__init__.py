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


VALID_HOMOSAPIENS_TR = get_json_resource("valid_homosapiens_tr.json")
HOMOSAPIENS_TR_SYNONYMS = get_json_resource("homosapiens_tr_synonyms.json")
HOMOSAPIENS_TR_AA_SEQUENCES = get_json_resource("homosapiens_tr_aa_sequences.json")
VALID_HOMOSAPIENS_MH = get_json_resource("valid_homosapiens_mh.json")
HOMOSAPIENS_MH_SYNONYMS = get_json_resource("homosapiens_mh_synonyms.json")


VALID_MUSMUSCULUS_TR = get_json_resource("valid_musmusculus_tr.json")
MUSMUSCULUS_TR_AA_SEQUENCES = get_json_resource("musmusculus_tr_aa_sequences.json")
VALID_MUSMUSCULUS_MH = get_json_resource("valid_musmusculus_mh.json")
MUSMUSCULUS_MH_SYNONYMS = get_json_resource("musmusculus_mh_synonyms.json")


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
