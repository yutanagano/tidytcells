from tidytcells._standardized_junction.standardized_junction import (
    StandardizedJunction,
)
from tidytcells._resources import MUSMUSCULUS_IG_AA_SEQUENCES


class StandardizedMusMusculusIgJunction(StandardizedJunction):
    _species = "musmusculus"
    _sequence_dictionary = MUSMUSCULUS_IG_AA_SEQUENCES
