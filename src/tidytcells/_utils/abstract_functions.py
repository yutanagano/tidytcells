from .gene_standardisers import GeneStandardiser
import re
from .._resources import AMINO_ACIDS
from typing import Dict, FrozenSet, Optional
from .warnings import *


def standardise_template(
    gene: str,
    gene_type: str,
    species: str,
    enforce_functional: bool,
    precision: str,
    suppress_warnings: bool,
    standardiser_dict: Dict[str, GeneStandardiser],
    allowed_precision: set,
) -> str:
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

    if not precision in allowed_precision:
        raise ValueError(f"precision must be in {allowed_precision}, got {precision}.")

    # For backward compatibility, fix CamelCased species
    species = "".join(species.split()).lower()

    if not species in standardiser_dict:
        if not suppress_warnings:
            warn_unsupported_species(species, gene_type)
        return gene

    original_input = gene

    gene = "".join(gene.split())
    gene = gene.replace("&nbsp;", "")
    gene = gene.upper()

    standardised = standardiser_dict[species](gene)

    invalid_reason = standardised.invalid(enforce_functional)
    if invalid_reason:
        if not suppress_warnings:
            warn_failure(
                reason_for_failure=invalid_reason,
                original_input=original_input,
                attempted_fix=standardised.compile("allele"),
                species=species,
            )
        return None

    return standardised.compile(precision)


def query_template(
    species: str,
    precision: str,
    functionality: str,
    contains: Optional[str],
    query_engine_dict: dict,
) -> FrozenSet[str]:
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
    if type(seq) != str:
        raise TypeError(f"seq must be type str, got {seq} ({type(seq)}).")

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
