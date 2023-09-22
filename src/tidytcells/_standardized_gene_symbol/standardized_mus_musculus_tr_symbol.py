from tidytcells._standardized_gene_symbol.standardized_tr_symbol import (
    StandardizedTrSymbol,
)
from tidytcells._resources import VALID_MUSMUSCULUS_TR


class StandardizedMusMusculusTrSymbol(StandardizedTrSymbol):
    _synonym_dictionary = dict()
    _valid_tr_dictionary = VALID_MUSMUSCULUS_TR
