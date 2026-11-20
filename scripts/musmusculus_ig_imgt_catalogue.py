import script_utility


def main() -> None:
    print("Fetching IG allele data from IMGT...")
    valid_allele_data = script_utility.get_ig_alleles_list("Mus+musculus")

    print("Fetching IG gene sequence data from IMGT...")
    sequence_data = script_utility.get_ig_aa_sequence_data("Mus+musculus")

    combined = script_utility.combine_symbols(valid_allele_data, sequence_data)
    json_output = {"provenance": script_utility.get_genedb_provenance(), "symbols": combined}

    script_utility.save_as_json(json_output, "musmusculus_ig_imgt_symbols.json")


if __name__ == "__main__":
    main()
