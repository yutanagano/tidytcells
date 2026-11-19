from tidytcells._query_engine.tr_query_engine import TrIgQueryEngine
from tidytcells._resources import HOMOSAPIENS_TR_IMGT_ALLELES


class HomoSapiensTrQueryEngine(TrIgQueryEngine):
    _valid_gene_dictionary = HOMOSAPIENS_TR_IMGT_ALLELES
