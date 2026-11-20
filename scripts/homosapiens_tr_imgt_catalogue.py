import script_utility


def main() -> None:
    print("Fetching TR allele data from IMGT...")
    valid_allele_data = script_utility.get_tr_alleles_list("Homo+sapiens")

    print("Fetching TR symbol synonyms from HGNC...")
    synonyms_data = script_utility.get_tr_synonyms_data(valid_allele_data)
    script_utility.save_as_json(synonyms_data, "homosapiens_tr_synonyms.json")

    print("Fetching TR gene sequence data from IMGT...")
    sequence_data = script_utility.get_tr_aa_sequence_data("Homo+sapiens")

    combined = script_utility.combine_symbols(valid_allele_data, sequence_data)
    json_output = {"provenance": script_utility.get_genedb_provenance(), "symbols": combined}

    script_utility.save_as_json(json_output, "homosapiens_tr_imgt_symbols.json")


if __name__ == "__main__":
    main()