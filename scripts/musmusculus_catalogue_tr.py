import script_utility


def main() -> None:
    print("Fetching TR allele data from IMGT...")
    valid_allele_data = script_utility.get_tr_alleles_list("Mus+musculus")
    script_utility.save_as_json(valid_allele_data, "valid_musmusculus_tr.json")


if __name__ == "__main__":
    main()