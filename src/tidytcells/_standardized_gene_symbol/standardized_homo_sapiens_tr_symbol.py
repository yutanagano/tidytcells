from tidytcells._standardized_gene_symbol.standardized_tr_symbol import (
    StandardizedTrSymbol,
)
from tidytcells._resources import VALID_HOMOSAPIENS_TR, HOMOSAPIENS_TR_SYNONYMS


class StandardizedHomoSapiensTrSymbol(StandardizedTrSymbol):
    _synonym_dictionary = HOMOSAPIENS_TR_SYNONYMS
    _valid_tr_dictionary = VALID_HOMOSAPIENS_TR
