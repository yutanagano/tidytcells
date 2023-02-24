'''
Warning message templates.
'''


from warnings import warn


def warn_failure(
    original_input: str,
    attempted_fix: str,
    species: str
):
    warn(
        f'Failed to standardise: "{original_input}" for species {species}. '
        f'Attempted fix "{attempted_fix}" did not meet the standardised '
        'format requirements. Ignoring this gene name...'
    )


def warn_unsupported_species(
    species: str,
    gene_type: str
):
    warn(
        f'Unsupported species: "{species}". '
        f'Skipping {gene_type} gene standardisation procedure...'
    )