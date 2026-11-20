from tidytcells._query_engine.tr_query_engine import TrIgQueryEngine
from tidytcells._resources import TR_SYMBOLS_BY_SP_DB


class HomoSapiensTrQueryEngine(TrIgQueryEngine):
    _valid_gene_dictionary = TR_SYMBOLS_BY_SP_DB["homosapiens"]["IMGT"]["symbols"]
