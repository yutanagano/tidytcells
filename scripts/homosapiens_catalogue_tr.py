import pandas as pd
from typing import Iterable
import script_utility


def main() -> None:
    print("Fetching TR allele data from IMGT...")
    valid_allele_data = script_utility.get_tr_alleles_list("Homo+sapiens")
    script_utility.save_as_json(valid_allele_data, "valid_homosapiens_tr.json")

    print("Fetching TR symbol synonyms from HGNC...")
    synonyms_data = get_synonyms_data(valid_allele_data)
    script_utility.save_as_json(synonyms_data, "homosapiens_tr_synonyms.json")

    print("Fetching TR gene sequence data from IMGT...")
    sequence_data = script_utility.get_tr_aa_sequence_data("Homo+sapiens")
    sequence_data["TRAJ35*01"]["J-CYS"] = "C"
    script_utility.save_as_json(sequence_data, "homosapiens_tr_aa_sequences.json")


def get_synonyms_data(valid_alleles: Iterable[str], is_tr: bool = True) -> dict:
    hgnc = script_utility.fetch_hgnc_data()

    field_name = "T cell receptor gene" if is_tr else "immunoglobulin gene"
    genes = hgnc[hgnc["Locus type"].str.contains(field_name)].copy()

    genes["Approved symbol"] = genes["Approved symbol"].str.replace(
        r"(?<!TR)DV", "/DV", regex=True
    )
    genes["Approved symbol"] = genes["Approved symbol"].str.replace(
        r"OR", "/OR", regex=True
    )

    # Only keep genes whose 'approved symbols' are in our IMGT list
    genes = genes[genes["Approved symbol"].map(lambda x: x in valid_alleles)].copy()

    # Get TR genes with "alias symbols"
    genes_with_aliases = genes[genes["Alias symbols"].notna()][
        ["Approved symbol", "Alias symbols"]
    ]
    genes_with_aliases["Alias symbols"] = genes_with_aliases["Alias symbols"].map(
        lambda x: [element.strip() for element in x.split(",")]
    )
    genes_with_aliases.columns = ["Approved symbol", "Synonym"]
    genes_with_aliases = genes_with_aliases.explode("Synonym")

    # Get TR genes with "previous symbols" (deprecated symbols)
    genes_with_depnames = genes[genes["Previous symbols"].notna()][
        ["Approved symbol", "Previous symbols"]
    ]
    genes_with_depnames["Previous symbols"] = genes_with_depnames[
        "Previous symbols"
    ].map(lambda x: [element.strip() for element in x.split(",")])
    genes_with_depnames.columns = ["Approved symbol", "Synonym"]
    genes_with_depnames = genes_with_depnames.explode("Synonym")

    # Combine
    synonyms = pd.concat([genes_with_aliases, genes_with_depnames])

    # Remove redundant synonyms
    synonyms = synonyms[synonyms["Approved symbol"] != synonyms["Synonym"]]

    # Remove ambiguous synonyms
    synonyms = synonyms.groupby("Synonym").aggregate(lambda x: x.tolist())
    synonyms = synonyms[synonyms["Approved symbol"].map(len) == 1].copy()
    synonyms["Approved symbol"] = synonyms["Approved symbol"].map(lambda x: x.pop())
    synonyms.index = synonyms.index.str.upper()

    # Remove any synonyms that are also names of other valid genes
    synonyms = synonyms[synonyms.index.map(lambda x: x not in valid_alleles)]

    synonyms_dict = synonyms["Approved symbol"].to_dict()

    # Remove ambiguous J1 and J2 (issue with both TCR and BCR)
    synonyms_dict.pop("J1", None)
    synonyms_dict.pop("J2", None)

    return synonyms_dict


if __name__ == "__main__":
    main()
