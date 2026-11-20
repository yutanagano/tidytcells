import json
from importlib.resources import files


def load_json(filename: str) -> dict:
    return json.loads(files(__name__).joinpath(filename).read_text())


TR_SYMBOLS_BY_SP_DB = {
    "homosapiens": {
        "IMGT": {
            "symbols": load_json("homosapiens_tr_imgt_symbols.json"),
            "synonyms": load_json("homosapiens_tr_synonyms.json"),
        },
    },
    "musmusculus": {
        "IMGT": {
            "symbols": load_json("musmusculus_tr_imgt_symbols.json"),
            "synonyms": dict(),
        },
    },
}


IG_SYMBOLS_BY_SP_DB = {
    "homosapiens": {
        "IMGT": {
            "symbols": load_json("homosapiens_ig_imgt_symbols.json"),
            "synonyms": load_json("homosapiens_ig_synonyms.json"),
        },
        "OGRDB": {
            "symbols": load_json("homosapiens_ig_ogrdb_symbols.json"),
            "synonyms": load_json("homosapiens_ig_synonyms.json"),
        },
    },
    "musmusculus": {
        "IMGT": {
            "symbols": load_json("musmusculus_ig_imgt_symbols.json"),
            "synonyms": dict(),
        },
    },
}


MH_SYMBOLS_BY_SP_DB = {
    "homosapiens": {
        "IPD": {
            "symbols": load_json("homosapiens_mh_ipd_symbols.json"),
            "synonyms_allele": dict(),
            "synonyms_gene": load_json("homosapiens_mh_ipd_synonyms_gene.json"),
        },
        "MRO": {
            "symbols": load_json("homosapiens_mh_mro_symbols.json"),
            "synonyms_allele": load_json("homosapiens_mh_mro_synonyms_allele.json"),
            "synonyms_gene": load_json("homosapiens_mh_mro_synonyms_gene.json"),
        },
    },
    "musmusculus": {
        "MRO": {
            "symbols": load_json("musmusculus_mh_mro_symbols.json"),
            "synonyms_allele": load_json("musmusculus_mh_mro_synonyms_allele.json"),
        },
    },
    "bostaurus": {
        "MRO": {
            "symbols": load_json("bostaurus_mh_mro_symbols.json"),
            "synonyms_allele": load_json("bostaurus_mh_mro_synonyms_allele.json"),
        },
    },
    "macacamulatta": {
        "MRO": {
            "symbols": load_json("macacamulatta_mh_mro_symbols.json"),
            "synonyms_allele": load_json("macacamulatta_mh_mro_synonyms_allele.json"),
        },
    },
    "pantroglodytes": {
        "MRO": {
            "symbols": load_json("pantroglodytes_mh_mro_symbols.json"),
            "synonyms_allele": load_json("pantroglodytes_mh_mro_synonyms_allele.json"),
        },
    },
    "rattusnorvegicus": {
        "MRO": {
            "symbols": load_json("rattusnorvegicus_mh_mro_symbols.json"),
            "synonyms_allele": load_json(
                "rattusnorvegicus_mh_mro_synonyms_allele.json"
            ),
        },
    },
    "susscrofa": {
        "MRO": {
            "symbols": load_json("susscrofa_mh_mro_symbols.json"),
            "synonyms_allele": load_json("susscrofa_mh_mro_synonyms_allele.json"),
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
