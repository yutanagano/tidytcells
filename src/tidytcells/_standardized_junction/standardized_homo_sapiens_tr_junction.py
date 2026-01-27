from tidytcells._standardized_junction.standardized_junction import (
    JunctionStandardizer,
)
from tidytcells._resources import HOMOSAPIENS_TR_AA_SEQUENCES


class HomoSapiensTrJunctionStandardizer(JunctionStandardizer):
    _species = "homosapiens"
    _sequence_dictionary = HOMOSAPIENS_TR_AA_SEQUENCES
