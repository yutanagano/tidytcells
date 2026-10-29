import json
from pathlib import Path

RESOURCES_DIR = Path(__file__).parent


def get_json_resource(filename: str) -> dict:
    with open(RESOURCES_DIR / filename, "r") as f:
        return json.load(f)


VALID_HOMOSAPIENS_TR = get_json_resource("valid_homosapiens_tr.json")
HOMOSAPIENS_TR_SYNONYMS = get_json_resource("homosapiens_tr_synonyms.json")
HOMOSAPIENS_TR_AA_SEQUENCES = get_json_resource("homosapiens_tr_aa_sequences.json")

VALID_HOMOSAPIENS_IG = get_json_resource("valid_homosapiens_ig.json")
HOMOSAPIENS_IG_SYNONYMS = get_json_resource("homosapiens_ig_synonyms.json")
HOMOSAPIENS_IG_AA_SEQUENCES = get_json_resource("homosapiens_ig_aa_sequences.json")

VALID_MUSMUSCULUS_TR = get_json_resource("valid_musmusculus_tr.json")
MUSMUSCULUS_TR_AA_SEQUENCES = get_json_resource("musmusculus_tr_aa_sequences.json")

VALID_MUSMUSCULUS_IG = get_json_resource("valid_musmusculus_ig.json")
MUSMUSCULUS_IG_AA_SEQUENCES = get_json_resource("musmusculus_ig_aa_sequences.json")



HOMOSAPIENS_MH_VALID_IPD = get_json_resource("ipd/homosapiens_mh_valid_allele_ipd.json")
HOMOSAPIENS_MH_SYNONYMS_GENE_HGNC_IPD = get_json_resource("ipd/homosapiens_mh_synonyms_hgnc.json")

HOMOSAPIENS_MH_VALID_PROTEIN_MRO = get_json_resource("mro/homosapiens_mh_valid_protein_mro.json")
HOMOSAPIENS_MH_SYNONYMS_ALLELE_MRO = get_json_resource("mro/homosapiens_mh_synonyms_mro.json")
HOMOSAPIENS_MH_SYNONYMS_GENE_HGNC_MRO = get_json_resource("mro/homosapiens_mh_synonyms_hgnc.json")

MUSMUSCULUS_MH_VALID_PROTEIN_MRO = get_json_resource("mro/musmusculus_mh_valid_protein_mro.json")
MUSMUSCULUS_MH_SYNONYMS_ALLELE_MRO = get_json_resource("mro/musmusculus_mh_synonyms_mro.json")

BOSTAURUS_MH_VALID_PROTEIN_MRO = get_json_resource("mro/bostaurus_mh_valid_protein_mro.json")
BOSTAURUS_MH_SYNONYMS_ALLELE_MRO = get_json_resource("mro/bostaurus_mh_synonyms_mro.json")

MACACAMULATTA_MH_VALID_PROTEIN_MRO = get_json_resource("mro/macacamulatta_mh_valid_protein_mro.json")
MACACAMULATTA_MH_SYNONYMS_ALLELE_MRO = get_json_resource("mro/macacamulatta_mh_synonyms_mro.json")

PANTROGLODYTES_MH_VALID_PROTEIN_MRO = get_json_resource("mro/pantroglodytes_mh_valid_protein_mro.json")
PANTROGLODYTES_MH_SYNONYMS_ALLELE_MRO = get_json_resource("mro/pantroglodytes_mh_synonyms_mro.json")

RATTUSNORVEGICUS_MH_VALID_PROTEIN_MRO = get_json_resource("mro/rattusnorvegicus_mh_valid_protein_mro.json")
RATTUSNORVEGICUS_MH_SYNONYMS_ALLELE_MRO = get_json_resource("mro/rattusnorvegicus_mh_synonyms_mro.json")

SUSSCROFA_MH_VALID_PROTEIN_MRO = get_json_resource("mro/susscrofa_mh_valid_protein_mro.json")
SUSSCROFA_MH_SYNONYMS_ALLELE_MRO = get_json_resource("mro/susscrofa_mh_synonyms_mro.json")


SUPPORTED_SPECIES_AND_THEIR_TR_AA_SEQUENCES = {
    "homosapiens": HOMOSAPIENS_TR_AA_SEQUENCES,
    "musmusculus": MUSMUSCULUS_TR_AA_SEQUENCES,
}

SUPPORTED_SPECIES_AND_THEIR_IG_AA_SEQUENCES = {
    "homosapiens": HOMOSAPIENS_IG_AA_SEQUENCES,
    "musmusculus": MUSMUSCULUS_IG_AA_SEQUENCES,
}

SUPPORTED_RECEPTOR_SPECIES_AND_THEIR_AA_SEQUENCES = {
    "TR": SUPPORTED_SPECIES_AND_THEIR_TR_AA_SEQUENCES,
    "IG": SUPPORTED_SPECIES_AND_THEIR_IG_AA_SEQUENCES
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
