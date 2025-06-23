from tidytcells._query_engine.tr_query_engine import TrIgQueryEngine
from tidytcells._resources import VALID_HOMOSAPIENS_TR


class HomoSapiensTrQueryEngine(TrIgQueryEngine):
    _valid_gene_dictionary = VALID_HOMOSAPIENS_TR
