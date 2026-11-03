import json
from pathlib import Path

RESOURCES_DIR = Path(__file__).parent


def get_json_resource(filename: str) -> dict:
    with open(RESOURCES_DIR / filename, "r") as f:
        return json.load(f)


VALID_HOMOSAPIENS_TR_IMGT = get_json_resource("valid_homosapiens_tr.json")
HOMOSAPIENS_TR_SYNONYMS = get_json_resource("homosapiens_tr_synonyms.json")
HOMOSAPIENS_TR_AA_SEQUENCES_IMGT = get_json_resource("homosapiens_tr_aa_sequences.json")
VALID_HOMOSAPIENS_MH = get_json_resource("valid_homosapiens_mh.json")
HOMOSAPIENS_MH_SYNONYMS = get_json_resource("homosapiens_mh_synonyms.json")
VALID_HOMOSAPIENS_MH_MRO = get_json_resource("valid_homosapiens_mh_mro.json")
HOMOSAPIENS_MH_SYNONYMS_ALLELE_MRO = get_json_resource("homosapiens_mh_synonyms_allele_mro.json")

VALID_HOMOSAPIENS_IG_IMGT = get_json_resource("valid_homosapiens_ig.json")
VALID_HOMOSAPIENS_IG_OGRDB = get_json_resource("valid_homosapiens_ig_ogrdb.json")
HOMOSAPIENS_IG_SYNONYMS = get_json_resource("homosapiens_ig_synonyms.json")
HOMOSAPIENS_IG_AA_SEQUENCES_IMGT = get_json_resource("homosapiens_ig_aa_sequences.json")
HOMOSAPIENS_IG_AA_SEQUENCES_OGRDB = get_json_resource("homosapiens_ig_aa_sequences_ogrdb.json")

VALID_MUSMUSCULUS_TR_IMGT = get_json_resource("valid_musmusculus_tr.json")
MUSMUSCULUS_TR_AA_SEQUENCES_IMGT = get_json_resource("musmusculus_tr_aa_sequences.json")
VALID_MUSMUSCULUS_MH = get_json_resource("valid_musmusculus_mh.json")
MUSMUSCULUS_MH_SYNONYMS = get_json_resource("musmusculus_mh_synonyms.json")
VALID_MUSMUSCULUS_MH_MRO = get_json_resource("valid_musmusculus_mh_mro.json")
MUSMUSCULUS_MH_SYNONYMS_ALLELE_MRO = get_json_resource("musmusculus_mh_synonyms_allele_mro.json")

VALID_MUSMUSCULUS_IG_IMGT = get_json_resource("valid_musmusculus_ig.json")
MUSMUSCULUS_IG_AA_SEQUENCES_IMGT = get_json_resource("musmusculus_ig_aa_sequences.json")




MRO_SPECIES_AND_VALID_SYMBOLS = {
    "homosapiens": HOMOSAPIENS_MH_VALID_PROTEIN_MRO,
    "bostaurus": BOSTAURUS_MH_VALID_PROTEIN_MRO,
    "macacamulatta": MACACAMULATTA_MH_VALID_PROTEIN_MRO,
    "pantroglodytes": PANTROGLODYTES_MH_VALID_PROTEIN_MRO,
    "susscrofa": SUSSCROFA_MH_VALID_PROTEIN_MRO,
    "musmusculus": MUSMUSCULUS_MH_VALID_MRO,
}

SUPPORTED_RECEPTOR_SPECIES_AND_THEIR_AA_SEQUENCES = {
    "homosapiens": {
        "TR": {"IMGT": HOMOSAPIENS_TR_AA_SEQUENCES_IMGT},
        "IG": {"IMGT": HOMOSAPIENS_IG_AA_SEQUENCES_IMGT, "OGRDB": HOMOSAPIENS_IG_AA_SEQUENCES_OGRDB},
    },
    "musmusculus": {
        "TR": {"IMGT": MUSMUSCULUS_TR_AA_SEQUENCES_IMGT},
        "IG": {"IMGT": MUSMUSCULUS_IG_AA_SEQUENCES_IMGT},
    }
}


SUPPORTED_IG_SPECIES_AND_THEIR_GENES = {
    "homosapiens": {
        "IMGT": {"genes": VALID_HOMOSAPIENS_IG_IMGT, "synonyms": HOMOSAPIENS_IG_SYNONYMS},
        "OGRDB": {"genes": VALID_HOMOSAPIENS_IG_OGRDB, "synonyms": HOMOSAPIENS_IG_SYNONYMS},
    },
    "musmusculus": {
        "IMGT": {"genes": VALID_MUSMUSCULUS_IG_IMGT, "synonyms": dict()},
    }
}

SUPPORTED_TR_SPECIES_AND_THEIR_GENES = {
    "homosapiens": {
        "IMGT": {"genes": VALID_HOMOSAPIENS_TR_IMGT, "synonyms": HOMOSAPIENS_TR_SYNONYMS},
    },
    "musmusculus": {
        "IMGT": {"genes": VALID_MUSMUSCULUS_TR_IMGT, "synonyms": dict()},
    }
}



HOMOSAPIENS_MH_VALID_IPD = get_json_resource("ipd_homosapiens_mh_valid_allele.json")
HOMOSAPIENS_MH_VALID_PROTEIN_MRO = get_json_resource("mro_homosapiens_mh_valid_protein.json")

MUSMUSCULUS_MH_VALID_PROTEIN_MRO = get_json_resource("mro_musmusculus_mh_valid_protein.json")

SUPPORTED_MH_SPECIES_AND_THEIR_GENES = {
    "homosapiens": {
        "IPD": {
            "valid_symbols": HOMOSAPIENS_MH_VALID_IPD,
            "allele_synonyms": dict(),
            "hgnc_synonyms": get_json_resource("ipd_homosapiens_mh_synonyms_hgnc.json")
        },
        "MRO": {
            "valid_symbols": HOMOSAPIENS_MH_VALID_PROTEIN_MRO,
            "allele_synonyms": get_json_resource("mro_homosapiens_mh_synonyms.json"),
            "hgnc_synonyms": get_json_resource("mro_homosapiens_mh_synonyms_hgnc.json")
        },
    },
    "musmusculus": {
        "MRO": {
            "valid_symbols": MUSMUSCULUS_MH_VALID_PROTEIN_MRO,
            "allele_synonyms": get_json_resource("mro_musmusculus_mh_synonyms.json")
        },
    },
    "bostaurus": {
        "MRO": {
            "valid_symbols": get_json_resource("mro_bostaurus_mh_valid_protein.json"),
            "allele_synonyms": get_json_resource("mro_bostaurus_mh_synonyms.json")
        },
    },
    "macacamulatta": {
        "MRO": {
            "valid_symbols": get_json_resource("mro_macacamulatta_mh_valid_protein.json"),
            "allele_synonyms": get_json_resource("mro_macacamulatta_mh_synonyms.json")
        },
    },
    "pantroglodytes": {
        "MRO": {
            "valid_symbols": get_json_resource("mro_pantroglodytes_mh_valid_protein.json"),
            "allele_synonyms": get_json_resource("mro_pantroglodytes_mh_synonyms.json")
        },
    },
    "rattusnorvegicus": {
        "MRO": {
            "valid_symbols": get_json_resource("mro_rattusnorvegicus_mh_valid_protein.json"),
            "allele_synonyms": get_json_resource("mro_rattusnorvegicus_mh_synonyms.json")
        },
    },
    "susscrofa": {
        "MRO": {
            "valid_symbols": get_json_resource("mro_susscrofa_mh_valid_protein.json"),
            "allele_synonyms": get_json_resource("mro_susscrofa_mh_synonyms.json")
        },
    },
}





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

AMINO_ACIDS_WITH_SPECIAL = AMINO_ACIDS | {"X", "*", "#"}
