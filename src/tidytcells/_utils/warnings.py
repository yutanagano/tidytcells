from logging import Logger


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


def warn_unsupported_species(species: str, gene_type: str, logger: Logger):
    logger.warning(
        f'Unsupported species: "{species}". ' f"Skipping {gene_type} standardisation."
    )
