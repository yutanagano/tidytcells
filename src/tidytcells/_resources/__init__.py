import json
from importlib.resources import files


def load_json(filename: str) -> dict:
    return json.loads(files(__name__).joinpath(filename).read_text())


HOMOSAPIENS_TR_IMGT_ALLELES = load_json("homosapiens_tr_imgt_alleles.json")
HOMOSAPIENS_TR_SYNONYMS = load_json("homosapiens_tr_synonyms.json")

HOMOSAPIENS_IG_IMGT_ALLELES = load_json("homosapiens_ig_imgt_alleles.json")
HOMOSAPIENS_IG_OGRDB_ALLELES = load_json("homosapiens_ig_ogrdb_alleles.json")
HOMOSAPIENS_IG_SYNONYMS = load_json("homosapiens_ig_synonyms.json")

MUSMUSCULUS_TR_IMGT_ALLELES = load_json("musmusculus_tr_imgt_alleles.json")

MUSMUSCULUS_IG_IMGT_ALLELES = load_json("musmusculus_ig_imgt_alleles.json")


SUPPORTED_RECEPTOR_SPECIES_AND_THEIR_GENES = {
    "homosapiens": {
        "TR": {
            "IMGT": {
                "genes": HOMOSAPIENS_TR_IMGT_ALLELES,
                "synonyms": HOMOSAPIENS_TR_SYNONYMS,
            },
        },
        "IG": {
            "IMGT": {
                "genes": HOMOSAPIENS_IG_IMGT_ALLELES,
                "synonyms": HOMOSAPIENS_IG_SYNONYMS,
            },
            "OGRDB": {
                "genes": HOMOSAPIENS_IG_OGRDB_ALLELES,
                "synonyms": HOMOSAPIENS_IG_SYNONYMS,
            },
        },
    },
    "musmusculus": {
        "TR": {
            "IMGT": {
                "genes": MUSMUSCULUS_TR_IMGT_ALLELES,
                "synonyms": dict(),
            },
        },
        "IG": {
            "IMGT": {
                "genes": MUSMUSCULUS_IG_IMGT_ALLELES,
                "synonyms": dict(),
            },
        },
    },
}


MH_SYMBOLS_BY_SP_DB = {
    "homosapiens": {
        "IPD": {
            "valid_symbols": load_json("homosapiens_mh_ipd_alleles.json"),
            "allele_synonyms": dict(),
            "hgnc_synonyms": load_json("homosapiens_mh_ipd_synonyms_gene.json"),
        },
        "MRO": {
            "valid_symbols": load_json("homosapiens_mh_mro_alleles.json"),
            "allele_synonyms": load_json("homosapiens_mh_mro_synonyms_allele.json"),
            "hgnc_synonyms": load_json("homosapiens_mh_mro_synonyms_gene.json"),
        },
    },
    "musmusculus": {
        "MRO": {
            "valid_symbols": load_json("musmusculus_mh_mro_alleles.json"),
            "allele_synonyms": load_json("musmusculus_mh_mro_synonyms_allele.json"),
        },
    },
    "bostaurus": {
        "MRO": {
            "valid_symbols": load_json("bostaurus_mh_mro_alleles.json"),
            "allele_synonyms": load_json("bostaurus_mh_mro_synonyms_allele.json"),
        },
    },
    "macacamulatta": {
        "MRO": {
            "valid_symbols": load_json("macacamulatta_mh_mro_alleles.json"),
            "allele_synonyms": load_json("macacamulatta_mh_mro_synonyms_allele.json"),
        },
    },
    "pantroglodytes": {
        "MRO": {
            "valid_symbols": load_json("pantroglodytes_mh_mro_alleles.json"),
            "allele_synonyms": load_json("pantroglodytes_mh_mro_synonyms_allele.json"),
        },
    },
    "rattusnorvegicus": {
        "MRO": {
            "valid_symbols": load_json("rattusnorvegicus_mh_mro_alleles.json"),
            "allele_synonyms": load_json(
                "rattusnorvegicus_mh_mro_synonyms_allele.json"
            ),
        },
    },
    "susscrofa": {
        "MRO": {
            "valid_symbols": load_json("susscrofa_mh_mro_alleles.json"),
            "allele_synonyms": load_json("susscrofa_mh_mro_synonyms_allele.json"),
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
