from tidytcells._standardized_gene_symbol.standardized_tr_symbol import (
    TrSymbolStandardizer,
)
from tidytcells._resources import VALID_HOMOSAPIENS_TR_IMGT, HOMOSAPIENS_TR_SYNONYMS


class HomoSapiensImgtTrSymbolStandardizer(TrSymbolStandardizer):
    _species = "homosapiens"
    _germline_database = "IMGT"
    _synonym_dictionary = {key: value for key, value in HOMOSAPIENS_TR_SYNONYMS.items() if value in VALID_HOMOSAPIENS_TR_IMGT}
    _valid_gene_dictionary = VALID_HOMOSAPIENS_TR_IMGT
    _valid_subgroups = {
        key.split("-")[0] for key in VALID_HOMOSAPIENS_TR_IMGT
    }