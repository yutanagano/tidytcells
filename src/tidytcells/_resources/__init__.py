import json
from importlib.resources import files


def load_json(filename: str) -> dict:
    return json.loads(files(__name__).joinpath(filename).read_text())


VALID_HOMOSAPIENS_TR_IMGT = load_json("homosapiens_tr_imgt_valid_genes.json")
HOMOSAPIENS_TR_SYNONYMS = load_json("homosapiens_tr_synonyms.json")
HOMOSAPIENS_TR_AA_SEQUENCES_IMGT = load_json("homosapiens_tr_imgt_aa_sequences.json")

VALID_HOMOSAPIENS_IG_IMGT = load_json("homosapiens_ig_imgt_valid_genes.json")
VALID_HOMOSAPIENS_IG_OGRDB = load_json("homosapiens_ig_ogrdb_valid_genes.json")
HOMOSAPIENS_IG_SYNONYMS = load_json("homosapiens_ig_synonyms.json")
HOMOSAPIENS_IG_AA_SEQUENCES_IMGT = load_json("homosapiens_ig_imgt_aa_sequences.json")
HOMOSAPIENS_IG_AA_SEQUENCES_OGRDB = load_json("homosapiens_ig_ogrdb_aa_sequences.json")

VALID_MUSMUSCULUS_TR_IMGT = load_json("musmusculus_tr_imgt_valid_genes.json")
MUSMUSCULUS_TR_AA_SEQUENCES_IMGT = load_json("musmusculus_tr_imgt_aa_sequences.json")

VALID_MUSMUSCULUS_IG_IMGT = load_json("musmusculus_ig_imgt_valid_genes.json")
MUSMUSCULUS_IG_AA_SEQUENCES_IMGT = load_json("musmusculus_ig_imgt_aa_sequences.json")


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
            "genes": VALID_HOMOSAPIENS_IG_IMGT,
            "synonyms": HOMOSAPIENS_IG_SYNONYMS,
        },
        "OGRDB": {
            "genes": VALID_HOMOSAPIENS_IG_OGRDB,
            "synonyms": HOMOSAPIENS_IG_SYNONYMS,
        },
    },
    "musmusculus": {
        "IMGT": {"genes": VALID_MUSMUSCULUS_IG_IMGT, "synonyms": dict()},
    },
}

SUPPORTED_TR_SPECIES_AND_THEIR_GENES = {
    "homosapiens": {
        "IMGT": {
            "genes": VALID_HOMOSAPIENS_TR_IMGT,
            "synonyms": HOMOSAPIENS_TR_SYNONYMS,
        },
    },
    "musmusculus": {
        "IMGT": {"genes": VALID_MUSMUSCULUS_TR_IMGT, "synonyms": dict()},
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
            "allele_synonyms": load_json("rattusnorvegicus_mh_mro_synonyms_allele.json"),
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
