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
    script_utility.save_as_json(sequence_data, "homosapiens_tr_aa_sequences.json")


def get_synonyms_data(valid_alleles: Iterable[str]) -> dict:
    hgnc = script_utility.fetch_hgnc_data()

    tr_genes = hgnc[hgnc["Locus type"].str.contains("T cell receptor gene")].copy()
    tr_genes["Approved symbol"] = tr_genes["Approved symbol"].str.replace(
        r"(?<!TR)DV", "/DV", regex=True
    )
    tr_genes["Approved symbol"] = tr_genes["Approved symbol"].str.replace(
        r"OR", "/OR", regex=True
    )

    # Only keep genes whose 'approved symbols' are in our IMGT list
    tr_genes = tr_genes[
        tr_genes["Approved symbol"].map(lambda x: x in valid_alleles)
    ].copy()

    # Get TR genes with "alias symbols"
    tr_genes_with_aliases = tr_genes[tr_genes["Alias symbols"].notna()][
        ["Approved symbol", "Alias symbols"]
    ]
    tr_genes_with_aliases["Alias symbols"] = tr_genes_with_aliases["Alias symbols"].map(
        lambda x: [element.strip() for element in x.split(",")]
    )
    tr_genes_with_aliases.columns = ["Approved symbol", "Synonym"]
    tr_genes_with_aliases = tr_genes_with_aliases.explode("Synonym")

    # Get TR genes with "previous symbols" (deprecated symbols)
    tr_genes_with_depnames = tr_genes[tr_genes["Previous symbols"].notna()][
        ["Approved symbol", "Previous symbols"]
    ]
    tr_genes_with_depnames["Previous symbols"] = tr_genes_with_depnames[
        "Previous symbols"
    ].map(lambda x: [element.strip() for element in x.split(",")])
    tr_genes_with_depnames.columns = ["Approved symbol", "Synonym"]
    tr_genes_with_depnames = tr_genes_with_depnames.explode("Synonym")

    # Combine
    tr_synonyms = pd.concat([tr_genes_with_aliases, tr_genes_with_depnames])

    # Remove redundant synonyms
    tr_synonyms = tr_synonyms[tr_synonyms["Approved symbol"] != tr_synonyms["Synonym"]]

    # Remove ambiguous synonyms
    tr_synonyms = tr_synonyms.groupby("Synonym").aggregate(lambda x: x.tolist())
    tr_synonyms = tr_synonyms[tr_synonyms["Approved symbol"].map(len) == 1].copy()
    tr_synonyms["Approved symbol"] = tr_synonyms["Approved symbol"].map(
        lambda x: x.pop()
    )
    tr_synonyms.index = tr_synonyms.index.str.upper()

    # Remove any synonyms that are also names of other valid genes
    tr_synonyms = tr_synonyms[tr_synonyms.index.map(lambda x: x not in valid_alleles)]

    return tr_synonyms["Approved symbol"].to_dict()


if __name__ == "__main__":
    main()
