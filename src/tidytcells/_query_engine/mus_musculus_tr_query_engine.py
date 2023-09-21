from tidytcells._query_engine.tr_query_engine import TrQueryEngine
from tidytcells._resources import MUSMUSCULUS_TCR


class MusMusculusTrQueryEngine(TrQueryEngine):
    _valid_tr_dictionary = MUSMUSCULUS_TCR
