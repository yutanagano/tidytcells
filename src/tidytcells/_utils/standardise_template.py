from .warnings import *


def standardise_template(
    gene: str,
    gene_type: str,
    species: str,
    enforce_functional: bool,
    precision: str,
    standardiser_dict: dict
) -> str:
    # Type errors
    if type(gene) != str:
        raise TypeError(
            'gene_name must be type str, got '
            f'{gene} ({type(gene)}).'
        )
    if type(species) != str:
        raise TypeError(
            'species must be type str, got '
            f'{species} ({type(species)}).'
        )
    if type(enforce_functional) != bool:
        raise TypeError(
            'enforce_functional must be type bool, got '
            f'{enforce_functional} ({type(enforce_functional)}).'
        )
    if type(precision) != str:
        raise TypeError(
            'precision must be type str, got '
            f'{precision} ({type(precision)}).'
        )
    
    # For backward compatibility, fix CamelCased species
    species = ''.join(species.split()).lower()

    # If the specified species is not supported, no-op (with warning)
    if not species in standardiser_dict:
        warn_unsupported_species(species, gene_type)
        return gene

    # Take note of initial input for reference
    original_input = gene

    # Clean whitespace, remove known pollutors, uppercase
    gene = ''.join(gene.split())
    gene = gene.replace('&nbsp;','')
    gene = gene.upper()

    # Standardisation attempt
    standardised = standardiser_dict[species](gene)

    if not standardised.valid(enforce_functional): # Standaridsation failure
        warn_failure(original_input, standardised.compile('allele'), species)
        return None
    
    return standardised.compile(precision)