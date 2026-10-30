from tidytcells._standardized_gene_symbol.standardized_ig_symbol import (
    IgSymbolStandardizer,
)
from tidytcells._resources import VALID_MUSMUSCULUS_IG_IMGT


class MusMusculusImgtIgSymbolStandardizer(IgSymbolStandardizer):
    _species = "musmusculus"
    _germline_database = "IMGT"
    _synonym_dictionary = dict()
    _valid_gene_dictionary = VALID_MUSMUSCULUS_IG_IMGT
    _valid_subgroups = {
        key.split("-")[0] for key in VALID_MUSMUSCULUS_IG_IMGT
    }