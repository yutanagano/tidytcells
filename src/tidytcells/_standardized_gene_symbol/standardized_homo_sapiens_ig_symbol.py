from tidytcells._standardized_gene_symbol.standardized_ig_symbol import (
    IgSymbolStandardizer,
)
from tidytcells._resources import VALID_HOMOSAPIENS_IG, HOMOSAPIENS_IG_SYNONYMS


class HomoSapiensIgSymbolStandardizer(IgSymbolStandardizer):
    _species = "homosapiens"
    _synonym_dictionary = HOMOSAPIENS_IG_SYNONYMS
    _valid_gene_dictionary = VALID_HOMOSAPIENS_IG
    _valid_subgroups = {
        key.split("-")[0] for key in VALID_HOMOSAPIENS_IG
    }