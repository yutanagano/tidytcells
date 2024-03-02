import script_utility


def main() -> None:
    print("Fetching TR allele data from IMGT...")
    valid_allele_data = script_utility.get_tr_alleles_list("Mus+musculus")
    script_utility.save_as_json(valid_allele_data, "valid_musmusculus_tr.json")

    print("Fetching TR gene sequence data from IMGT...")
    sequence_data = script_utility.get_tr_aa_sequence_data("Mus+musculus")
    script_utility.save_as_json(sequence_data, "musmusculus_tr_aa_sequences.json")


if __name__ == "__main__":
    main()
