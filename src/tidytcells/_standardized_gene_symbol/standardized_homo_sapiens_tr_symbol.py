from tidytcells._standardized_gene_symbol.standardized_tr_symbol import (
    StandardizedTrSymbol,
)
from tidytcells._resources import HOMOSAPIENS_TCR, HOMOSAPIENS_TCR_SYNONYMS


class StandardizedHomoSapiensTrSymbol(StandardizedTrSymbol):
    _synonym_dictionary = HOMOSAPIENS_TCR_SYNONYMS
    _valid_tcr_dictionary = HOMOSAPIENS_TCR
