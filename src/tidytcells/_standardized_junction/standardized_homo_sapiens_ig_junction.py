from tidytcells._standardized_junction.standardized_junction import (
    JunctionStandardizer,
)
from tidytcells._resources import HOMOSAPIENS_IG_AA_SEQUENCES


class HomoSapiensIgJunctionStandardizer(JunctionStandardizer):
    _species = "homosapiens"
    _sequence_dictionary = HOMOSAPIENS_IG_AA_SEQUENCES
