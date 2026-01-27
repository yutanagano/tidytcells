from tidytcells._standardized_gene_symbol.standardized_ig_symbol import (
    IgSymbolStandardizer,
)
from tidytcells._resources import VALID_MUSMUSCULUS_IG


class MusMusculusIgSymbolStandardizer(IgSymbolStandardizer):
    _species = "musmusculus"
    _synonym_dictionary = dict()
    _valid_gene_dictionary = VALID_MUSMUSCULUS_IG
    _valid_subgroups = {
        key.split("-")[0] for key in VALID_MUSMUSCULUS_IG
    }