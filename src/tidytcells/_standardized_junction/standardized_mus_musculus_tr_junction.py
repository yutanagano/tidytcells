from tidytcells._standardized_junction.standardized_junction import (
    JunctionStandardizer,
)
from tidytcells._resources import MUSMUSCULUS_TR_AA_SEQUENCES


class MusMusculusTrJunctionStandardizer(JunctionStandardizer):
    _species = "musmusculus"
    _sequence_dictionary = MUSMUSCULUS_TR_AA_SEQUENCES
