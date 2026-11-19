from tidytcells._query_engine.tr_query_engine import TrIgQueryEngine
from tidytcells._resources import MUSMUSCULUS_TR_IMGT_ALLELES


class MusMusculusTrQueryEngine(TrIgQueryEngine):
    _valid_gene_dictionary = MUSMUSCULUS_TR_IMGT_ALLELES
