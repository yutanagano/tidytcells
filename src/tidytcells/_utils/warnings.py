from logging import Logger
from typing import Union

from tidytcells.result import MhGene, ReceptorGene, Junction


def warn_failure(
    reason_for_failure: str,
    original_input: str,
    attempted_fix: str,
    species: str,
    logger: Logger,
):
    warning_message = f'Failed to standardize "{original_input}" for species {species}: {reason_for_failure}.'

    if original_input != attempted_fix:
        warning_message += f' (best attempted fix: "{attempted_fix}").'

    logger.warning(warning_message)


def warn_result_failure(
    result: Union[MhGene, ReceptorGene, Junction],
    logger: Logger,
):
    warn_failure(reason_for_failure=result.error,
                 original_input=result.original_input,
                 attempted_fix=result.attempted_fix,
                 species=result.species,
                 logger=logger)

def warn_unsupported_species(species: str, type: str, logger: Logger):
    logger.warning(
        f'Unsupported species: "{species}". ' f"Skipping {type} standardization."
    )
