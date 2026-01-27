from tidytcells._standardized_gene_symbol.standardized_tr_symbol import (
    TrSymbolStandardizer,
)
from tidytcells._resources import VALID_HOMOSAPIENS_TR, HOMOSAPIENS_TR_SYNONYMS


class HomoSapiensTrSymbolStandardizer(TrSymbolStandardizer):
    _species = "homosapiens"
    _synonym_dictionary = HOMOSAPIENS_TR_SYNONYMS
    _valid_gene_dictionary = VALID_HOMOSAPIENS_TR
    _valid_subgroups = {
        key.split("-")[0] for key in VALID_HOMOSAPIENS_TR
    }