import script_utility


def main() -> None:
    print("Fetching TR allele data from IMGT...")
    valid_allele_data = script_utility.get_tr_alleles_list("Mus+musculus")

    print("Fetching TR gene sequence data from IMGT...")
    sequence_data = script_utility.get_tr_aa_sequence_data("Mus+musculus")

    combined = script_utility.combine_symbols(valid_allele_data, sequence_data)
    script_utility.save_as_json(combined, "musmusculus_tr_imgt_symbols.json")


if __name__ == "__main__":
    main()
