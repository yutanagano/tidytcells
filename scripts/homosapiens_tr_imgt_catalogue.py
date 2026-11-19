import script_utility


def main() -> None:
    print("Fetching TR allele data from IMGT...")
    valid_allele_data = script_utility.get_tr_alleles_list("Homo+sapiens")
    script_utility.save_as_json(
        valid_allele_data, "homosapiens_tr_imgt_valid_genes.json"
    )

    print("Fetching TR symbol synonyms from HGNC and legacy mAb names...")
    synonyms_data = {
        **script_utility.get_tr_hgnc_synonyms_data(valid_allele_data),
        **script_utility.get_tr_legacy_synonyms_data("Homo sapiens", valid_allele_data),
    }
    script_utility.save_as_json(synonyms_data, "homosapiens_tr_synonyms.json")

    print("Fetching TR gene sequence data from IMGT...")
    sequence_data = script_utility.get_tr_aa_sequence_data("Homo+sapiens")
    script_utility.save_as_json(sequence_data, "homosapiens_tr_imgt_aa_sequences.json")


if __name__ == "__main__":
    main()
