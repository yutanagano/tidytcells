from tidytcells._standardized_gene_symbol.standardized_ig_symbol import (
    IgSymbolStandardizer,
)
from tidytcells._resources import VALID_HOMOSAPIENS_IG_IMGT, VALID_HOMOSAPIENS_IG_OGRDB, HOMOSAPIENS_IG_SYNONYMS


class HomoSapiensImgtIgSymbolStandardizer(IgSymbolStandardizer):
    _species = "homosapiens"
    _germline_database = "IMGT"
    _synonym_dictionary = {key: value for key, value in HOMOSAPIENS_IG_SYNONYMS.items() if value in VALID_HOMOSAPIENS_IG_IMGT}
    _valid_gene_dictionary = VALID_HOMOSAPIENS_IG_IMGT
    _valid_subgroups = {
        key.split("-")[0] for key in VALID_HOMOSAPIENS_IG_IMGT
    }

class HomoSapiensOgrdbIgSymbolStandardizer(IgSymbolStandardizer):
    _species = "homosapiens"
    _germline_database = "OGRDB"
    _synonym_dictionary = {key: value for key, value in HOMOSAPIENS_IG_SYNONYMS.items() if value in VALID_HOMOSAPIENS_IG_OGRDB}
    _valid_gene_dictionary = VALID_HOMOSAPIENS_IG_OGRDB
    _valid_subgroups = {
        key.split("-")[0] for key in VALID_HOMOSAPIENS_IG_OGRDB
    }