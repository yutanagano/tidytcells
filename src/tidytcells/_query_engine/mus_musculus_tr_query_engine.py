from tidytcells._query_engine.tr_query_engine import TrIgQueryEngine
from tidytcells._resources import SUPPORTED_RECEPTOR_SPECIES_AND_THEIR_GENES


class MusMusculusTrQueryEngine(TrIgQueryEngine):
    _valid_gene_dictionary = SUPPORTED_RECEPTOR_SPECIES_AND_THEIR_GENES["musmusculus"][
        "TR"
    ]["IMGT"]["symbols"]
