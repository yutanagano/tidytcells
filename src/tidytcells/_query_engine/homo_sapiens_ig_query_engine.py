from tidytcells._query_engine.tr_query_engine import TrIgQueryEngine
from tidytcells._resources import HOMOSAPIENS_IG_IMGT_ALLELES


class HomoSapiensIgQueryEngine(TrIgQueryEngine):
    _valid_gene_dictionary = HOMOSAPIENS_IG_IMGT_ALLELES
