from tidytcells._query_engine.tr_query_engine import TrQueryEngine
from tidytcells._resources import VALID_HOMOSAPIENS_IG


class HomoSapiensIgQueryEngine(TrQueryEngine):
    _valid_tr_dictionary = VALID_HOMOSAPIENS_IG
