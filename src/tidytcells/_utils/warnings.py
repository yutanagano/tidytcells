from warnings import warn


def warn_failure(
    reason_for_failure: str, original_input: str, attempted_fix: str, species: str
):
    warning_message = f'Failed to standardize "{original_input}" for species {species}: {reason_for_failure}.'

    if original_input != attempted_fix:
        warning_message += f' (best attempted fix: "{attempted_fix}").'

    warn(warning_message)


def warn_unsupported_species(species: str, gene_type: str):
    warn(f'Unsupported species: "{species}". ' f"Skipping {gene_type} standardisation.")
