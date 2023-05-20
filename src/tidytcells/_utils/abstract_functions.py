import re
from .._resources import AMINO_ACIDS
from typing import FrozenSet, Optional
from .warnings import *


def standardise_template(
    gene: str,
    gene_type: str,
    species: str,
    enforce_functional: bool,
    precision: str,
    suppress_warnings: bool,
    standardiser_dict: dict,
    allowed_precision: set,
) -> str:
    # Type checks
    if type(gene) != str:
        raise TypeError(f"gene_name must be type str, got {gene} ({type(gene)}).")
    if type(species) != str:
        raise TypeError(f"species must be type str, got {species} ({type(species)}).")
    if type(enforce_functional) != bool:
        raise TypeError(
            "enforce_functional must be type bool, got "
            f"{enforce_functional} ({type(enforce_functional)})."
        )
    if type(precision) != str:
        raise TypeError(
            f"precision must be type str, got {precision} ({type(precision)})."
        )
    if type(suppress_warnings) != bool:
        raise TypeError(
            "suppress_warnings must be type bool, got "
            f"{suppress_warnings} ({type(suppress_warnings)})."
        )

    # Allowed value checks
    if not precision in allowed_precision:
        raise ValueError(f"precision must be in {allowed_precision}, got {precision}.")

    # For backward compatibility, fix CamelCased species
    species = "".join(species.split()).lower()

    # If the specified species is not supported, no-op (with warning)
    if not species in standardiser_dict:
        if not suppress_warnings:
            warn_unsupported_species(species, gene_type)
        return gene

    # Take note of initial input for reference
    original_input = gene

    # Clean whitespace, remove known pollutors, uppercase
    gene = "".join(gene.split())
    gene = gene.replace("&nbsp;", "")
    gene = gene.upper()

    # Standardisation attempt
    standardised = standardiser_dict[species](gene)

    if not standardised.valid(enforce_functional):  # Standaridsation failure
        if not suppress_warnings:
            warn_failure(original_input, standardised.compile("allele"), species)
        return None

    return standardised.compile(precision)


def query_template(
    species: str,
    precision: str,
    functionality: str,
    contains: Optional[str],
    query_engine_dict: dict,
) -> FrozenSet[str]:
    # Type checks
    if type(species) != str:
        raise TypeError(f"species must be type str, got {species} ({type(species)}).")
    if type(functionality) != str:
        raise TypeError(
            f"functionality must be type str, got {functionality} ({type(functionality)})."
        )
    if type(precision) != str:
        raise TypeError(
            f"precision must be type str, got {precision} ({type(precision)})."
        )
    if not (contains is None or type(contains) == str):
        raise TypeError(
            f"contains must be either None or type str, got {contains} ({type(contains)})."
        )

    # Allowed value checks
    if not precision in {"allele", "gene"}:
        raise ValueError(
            f'precision must be either "allele" or "gene", got {precision}.'
        )
    if not functionality in {"any", "F", "NF", "P", "ORF"}:
        raise ValueError(
            f'functionality must be "any", "F", "NF", "P", or "ORF", got {functionality}.'
        )

    # For backward compatibility, fix CamelCased species
    species = "".join(species.split()).lower()

    if not species in query_engine_dict:
        raise ValueError(f"Unsupported species: {species}. No data available.")

    result = query_engine_dict[species].query(
        precision=precision, functionality=functionality
    )

    if contains is None:
        return result

    return frozenset([i for i in result if re.search(contains, i)])


def standardise_aa_template(seq: str, suppress_warnings: bool):
    # type checks
    if type(seq) != str:
        raise TypeError(f"seq must be type str, got {seq} ({type(seq)}).")

    # take note of original input
    original_input = seq

    seq = seq.upper()

    for char in seq:
        if not char in AMINO_ACIDS:
            if not suppress_warnings:
                warn(
                    f"Input {original_input} was rejected as it is not a valid amino acid sequence."
                )
            return None

    return seq
