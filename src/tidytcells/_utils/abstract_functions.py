from .gene_standardizers import GeneStandardizer
import re
from .._resources import AMINO_ACIDS
from typing import Dict, FrozenSet, Optional
from .warnings import *


def standardize_template(
    gene: str,
    gene_type: str,
    species: str,
    enforce_functional: bool,
    precision: str,
    on_fail: str,
    suppress_warnings: bool,
    standardizer_dict: Dict[str, GeneStandardizer],
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
    if type(on_fail) != str:
        raise TypeError(f"on_fail must be type str, got {on_fail} ({type(on_fail)}).")
    if type(suppress_warnings) != bool:
        raise TypeError(
            "suppress_warnings must be type bool, got "
            f"{suppress_warnings} ({type(suppress_warnings)})."
        )

    if not precision in allowed_precision:
        raise ValueError(f"precision must be in {allowed_precision}, got {precision}.")
    if not on_fail in ("reject", "keep"):
        raise ValueError(f'on_fail must be "reject" or "keep", got {on_fail}.')

    # For backward compatibility, fix CamelCased species
    species = "".join(species.split()).lower()

    if not species in standardizer_dict:
        if not suppress_warnings:
            warn_unsupported_species(species, gene_type)
        return gene

    standardized = standardizer_dict[species](gene)

    invalid_reason = standardized.invalid(enforce_functional)
    if invalid_reason:
        if not suppress_warnings:
            warn_failure(
                reason_for_failure=invalid_reason,
                original_input=gene,
                attempted_fix=standardized.compile("allele"),
                species=species,
            )
        if on_fail == "reject":
            return None
        return gene

    return standardized.compile(precision)


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


def standardize_aa_template(seq: str, on_fail: str, suppress_warnings: bool):
    if type(seq) != str:
        raise TypeError(f"seq must be type str, got {seq} ({type(seq)}).")
    if type(on_fail) != str:
        raise TypeError(f"on_fail must be type str, got {on_fail} ({type(on_fail)}).")
    if type(suppress_warnings) != bool:
        raise TypeError(
            "suppress_warnings must be type bool, got "
            f"{suppress_warnings} ({type(suppress_warnings)})."
        )

    if not on_fail in ("reject", "keep"):
        raise ValueError(f'on_fail must be "reject" or "keep", got {on_fail}.')

    original_input = seq

    seq = seq.upper()

    for char in seq:
        if not char in AMINO_ACIDS:
            if not suppress_warnings:
                warn(
                    f"Failed to standardize {original_input}: not a valid amino acid sequence."
                )
            if on_fail == "reject":
                return None
            return original_input

    return seq
