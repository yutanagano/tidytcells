from tidytcells._query_engine.tr_query_engine import TrQueryEngine
from tidytcells._resources import HOMOSAPIENS_TCR


class HomoSapiensTrQueryEngine(TrQueryEngine):
    _valid_tr_dictionary = HOMOSAPIENS_TCR
