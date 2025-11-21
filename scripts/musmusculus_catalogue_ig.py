import script_utility
from homosapiens_catalogue_tr import get_synonyms_data


def main() -> None:
    print("Fetching IG allele data from IMGT...")
    valid_allele_data = script_utility.get_ig_alleles_list("Mus+musculus")
    script_utility.save_as_json(valid_allele_data, "valid_musmusculus_ig.json")

    print("Fetching IG gene sequence data from IMGT...")
    sequence_data = script_utility.get_ig_aa_sequence_data("Mus+musculus")
    script_utility.save_as_json(sequence_data, "musmusculus_ig_aa_sequences.json")


if __name__ == "__main__":
    main()
