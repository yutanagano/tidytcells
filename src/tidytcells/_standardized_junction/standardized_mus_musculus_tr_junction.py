from tidytcells._standardized_junction.standardized_junction import (
    StandardizedJunction,
)
from tidytcells._resources import MUSMUSCULUS_TR_AA_SEQUENCES


class StandardizedMusMusculusTrJunction(StandardizedJunction):
    _species = "musmusculus"
    _sequence_dictionary = MUSMUSCULUS_TR_AA_SEQUENCES
