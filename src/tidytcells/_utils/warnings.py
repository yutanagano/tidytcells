from logging import Logger
from typing import Optional, Union

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

def warn_unsupported_species(species: str, module: str, valid_species, logger: Logger):
    logger.warning(
        f'Unsupported species: "{species}" (valid options are: {sorted(valid_species)}). '
        f"Skipping {module} standardization."
    )


def warn_unsupported_receptor_database(
    receptor_database: str,
    species: str,
    module: str,
    valid_databases,
    logger: Logger,
    log_failures: bool = True,
    locus: Optional[str] = None,
) -> str:
    '''
    Log a warning for an unsupported receptor_database setting (if `log_failures`)
    and return the corresponding error message.
    '''
    locus_part = f' and locus "{locus}"' if locus is not None else ""
    error = (
        f'Unsupported receptor_database: "{receptor_database}" for species "{species}"{locus_part} '
        f'(valid options are: {sorted(valid_databases)}). '
        f'Use tidytcells.config.set_receptor_database() to switch.'
    )
    if log_failures:
        logger.warning(f'{error} Skipping {module} standardization.')

    return error
