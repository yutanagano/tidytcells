from tidytcells._standardized_junction.standardized_junction import (
    StandardizedJunction,
)
from tidytcells._resources import HOMOSAPIENS_IG_AA_SEQUENCES


class StandardizedHomoSapiensIgJunction(StandardizedJunction):
    _sequence_dictionary = HOMOSAPIENS_IG_AA_SEQUENCES
