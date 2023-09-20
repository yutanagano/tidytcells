from tidytcells._standardized_gene_symbol.standardized_tr_symbol import (
    StandardizedTrSymbol,
)
from tidytcells._resources import MUSMUSCULUS_TCR


class StandardizedMusMusculusTrSymbol(StandardizedTrSymbol):
    _synonym_dictionary = dict()
    _valid_tcr_dictionary = MUSMUSCULUS_TCR
