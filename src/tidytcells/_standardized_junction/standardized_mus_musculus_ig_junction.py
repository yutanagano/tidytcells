from tidytcells._standardized_junction.standardized_junction import (
    JunctionStandardizer,
)
from tidytcells._resources import MUSMUSCULUS_IG_AA_SEQUENCES


class MusMusculusIgJunctionStandardizer(JunctionStandardizer):
    _species = "musmusculus"
    _sequence_dictionary = MUSMUSCULUS_IG_AA_SEQUENCES
