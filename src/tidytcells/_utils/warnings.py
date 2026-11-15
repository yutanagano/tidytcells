from typing import Optional, Union

from tidytcells.result import MhGene, ReceptorGene, Junction


def unsupported_species_message(species: str, valid_species) -> str:
    """
    Format an error message for an unsupported species.
    """
    return (
        f'Unsupported species: "{species}" '
        f"(valid options are: {sorted(valid_species)})."
    )


def unsupported_receptor_database_message(
    receptor_database: str,
    species: str,
    valid_databases,
    locus: Optional[str] = None,
) -> str:
    """
    Format an error message for an unsupported receptor_database setting.
    """
    locus_part = f' and locus "{locus}"' if locus is not None else ""
    return (
        f'Unsupported receptor_database: "{receptor_database}" for species "{species}"{locus_part} '
        f"(valid options are: {sorted(valid_databases)}). "
        "Use tidytcells.config.set_receptor_database() to switch."
    )


def unsupported_mh_database_message(
    mh_database: str,
    species: str,
    valid_databases,
) -> str:
    """
    Format an error message for an unsupported mh_database setting.
    """
    return (
        f'Unsupported mh_database: "{mh_database}" for species "{species}" '
        f"(valid options are: {sorted(valid_databases)}). "
        "Use tidytcells.config.set_mh_database() to switch."
    )


def unsupported_locus_message(locus: str, species: str, valid_loci) -> str:
    """
    Format an error message for an unsupported locus.
    """
    return (
        f'Unsupported locus: "{locus}" for species "{species}" '
        f"(valid options are: {sorted(valid_loci)})."
    )


def standardization_failure_message(
    result: Union[MhGene, ReceptorGene, Junction],
) -> str:
    """
    Format a warning message for a failed standardization result.
    """
    message = (
        f'Failed to standardize "{result.original_input}" '
        f"for species {result.species}: {result.error}."
    )
    if result.original_input != result.attempted_fix:
        message += f' (best attempted fix: "{result.attempted_fix}").'
    return message
