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


def _flat_sequence_dict(allele_tree: dict) -> dict:
    return {
        f"{gene}*{allele}": dict(data)
        for gene, alleles in allele_tree.items()
        for allele, data in alleles.items()
        if any(field != "functionality" for field in data)
    }


HOMOSAPIENS_TR_AA_SEQUENCES_IMGT = _flat_sequence_dict(HOMOSAPIENS_TR_IMGT_ALLELES)

HOMOSAPIENS_IG_AA_SEQUENCES_IMGT = _flat_sequence_dict(HOMOSAPIENS_IG_IMGT_ALLELES)
HOMOSAPIENS_IG_AA_SEQUENCES_OGRDB = _flat_sequence_dict(HOMOSAPIENS_IG_OGRDB_ALLELES)

MUSMUSCULUS_TR_AA_SEQUENCES_IMGT = _flat_sequence_dict(MUSMUSCULUS_TR_IMGT_ALLELES)

MUSMUSCULUS_IG_AA_SEQUENCES_IMGT = _flat_sequence_dict(MUSMUSCULUS_IG_IMGT_ALLELES)


SUPPORTED_RECEPTOR_SPECIES_AND_THEIR_AA_SEQUENCES = {
    "homosapiens": {
        "TR": {"IMGT": HOMOSAPIENS_TR_AA_SEQUENCES_IMGT},
        "IG": {
            "IMGT": HOMOSAPIENS_IG_AA_SEQUENCES_IMGT,
            "OGRDB": HOMOSAPIENS_IG_AA_SEQUENCES_OGRDB,
        },
    },
    "musmusculus": {
        "TR": {"IMGT": MUSMUSCULUS_TR_AA_SEQUENCES_IMGT},
        "IG": {"IMGT": MUSMUSCULUS_IG_AA_SEQUENCES_IMGT},
    },
}


SUPPORTED_IG_SPECIES_AND_THEIR_GENES = {
    "homosapiens": {
        "IMGT": {
            "genes": HOMOSAPIENS_IG_IMGT_ALLELES,
            "synonyms": HOMOSAPIENS_IG_SYNONYMS,
        },
        "OGRDB": {
            "genes": HOMOSAPIENS_IG_OGRDB_ALLELES,
            "synonyms": HOMOSAPIENS_IG_SYNONYMS,
        },
    },
    "musmusculus": {
        "IMGT": {"genes": MUSMUSCULUS_IG_IMGT_ALLELES, "synonyms": dict()},
    },
}

SUPPORTED_TR_SPECIES_AND_THEIR_GENES = {
    "homosapiens": {
        "IMGT": {
            "genes": HOMOSAPIENS_TR_IMGT_ALLELES,
            "synonyms": HOMOSAPIENS_TR_SYNONYMS,
        },
    },
    "musmusculus": {
        "IMGT": {"genes": MUSMUSCULUS_TR_IMGT_ALLELES, "synonyms": dict()},
    },
}


MH_SYMBOLS_BY_SP_DB = {
    "homosapiens": {
        "IPD": {
            "valid_symbols": load_json("homosapiens_mh_ipd_valid_alleles.json"),
            "allele_synonyms": dict(),
            "hgnc_synonyms": load_json("homosapiens_mh_ipd_synonyms_gene.json"),
        },
        "MRO": {
            "valid_symbols": load_json("homosapiens_mh_mro_valid_proteins.json"),
            "allele_synonyms": load_json("homosapiens_mh_mro_synonyms_allele.json"),
            "hgnc_synonyms": load_json("homosapiens_mh_mro_synonyms_gene.json"),
        },
    },
    "musmusculus": {
        "MRO": {
            "valid_symbols": load_json("musmusculus_mh_mro_valid_proteins.json"),
            "allele_synonyms": load_json("musmusculus_mh_mro_synonyms_allele.json"),
        },
    },
    "bostaurus": {
        "MRO": {
            "valid_symbols": load_json("bostaurus_mh_mro_valid_proteins.json"),
            "allele_synonyms": load_json("bostaurus_mh_mro_synonyms_allele.json"),
        },
    },
    "macacamulatta": {
        "MRO": {
            "valid_symbols": load_json("macacamulatta_mh_mro_valid_proteins.json"),
            "allele_synonyms": load_json("macacamulatta_mh_mro_synonyms_allele.json"),
        },
    },
    "pantroglodytes": {
        "MRO": {
            "valid_symbols": load_json("pantroglodytes_mh_mro_valid_proteins.json"),
            "allele_synonyms": load_json("pantroglodytes_mh_mro_synonyms_allele.json"),
        },
    },
    "rattusnorvegicus": {
        "MRO": {
            "valid_symbols": load_json("rattusnorvegicus_mh_mro_valid_proteins.json"),
            "allele_synonyms": load_json(
                "rattusnorvegicus_mh_mro_synonyms_allele.json"
            ),
        },
    },
    "susscrofa": {
        "MRO": {
            "valid_symbols": load_json("susscrofa_mh_mro_valid_proteins.json"),
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
