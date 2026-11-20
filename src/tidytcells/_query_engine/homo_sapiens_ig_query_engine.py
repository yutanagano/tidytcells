from tidytcells._query_engine.tr_query_engine import TrIgQueryEngine
from tidytcells._resources import IG_SYMBOLS_BY_SP_DB


class HomoSapiensIgQueryEngine(TrIgQueryEngine):
    _valid_gene_dictionary = IG_SYMBOLS_BY_SP_DB["homosapiens"]["IMGT"]["symbols"]
