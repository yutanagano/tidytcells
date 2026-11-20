import script_utility


def main() -> None:
    print("Fetching IG allele data from IMGT...")
    valid_allele_data = script_utility.get_ig_alleles_list("Homo+sapiens")

    print("Fetching IG symbol synonyms from HGNC...")
    synonyms_data = script_utility.get_ig_synonyms_data(valid_allele_data)
    script_utility.save_as_json(synonyms_data, "homosapiens_ig_synonyms.json")

    print("Fetching IG gene sequence data from IMGT...")
    sequence_data = script_utility.get_ig_aa_sequence_data("Homo+sapiens")

    combined = script_utility.combine_symbols(valid_allele_data, sequence_data)
    json_output = {"provenance": script_utility.get_genedb_provenance(), "symbols": combined}

    script_utility.save_as_json(json_output, "homosapiens_ig_imgt_symbols.json")


if __name__ == "__main__":
    main()
