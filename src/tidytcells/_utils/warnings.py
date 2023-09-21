from warnings import warn


def warn_failure(
    reason_for_failure: str, original_input: str, attempted_fix: str, species: str
):
    warn(
        f'Failed to standardize "{original_input}" for species {species}: '
        f'{reason_for_failure}. Attempted fix: "{attempted_fix}".'
    )


def warn_unsupported_species(species: str, gene_type: str):
    warn(f'Unsupported species: "{species}". ' f"Skipping {gene_type} standardisation.")
