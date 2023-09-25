# Valid values for TR and MH genes

## Homo sapiens

### Valid TR genes

- Data source: IMGT/GENE-DB (https://www.imgt.org/genedb/)
    - Use interactive GENE-DB interface to retrieve all Homo sapiens TR genes
    - Then, in order to get a list of all known alleles of each gene, request the nucleotide sequence of all functional, ORF and pseudogene alleles in FASTA format
    - Downloaded FASTA data found at: https://github.com/yutanagano/tidytcells/tree/main/data/homosapiens_tr.fasta
- Last date of access: 22nd September 2023
- Data parsing:
    - Using the FASTA data retrieved as above, use the data processing script (found at: https://github.com/yutanagano/tidytcells/tree/main/scripts/homosapiens_catalogue_tr.py) to transform the data to json form, where the data is encoded as:
    ```
    {
        gene name: {
            allele number: functionality
        }
    }
    ```

### Valid MH genes

- Data source: IMGT/MH-DB (https://github.com/ANHIG/IMGTHLA/blob/Latest/wmda/)
    - Directly pull HLA nomenclature files from the web using github's raw file API (e.g. https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom.txt)
- Last date of access: 22nd September 2023
- Data parsing:
    - The data processing script (found at: https://github.com/yutanagano/tidytcells/tree/main/scripts/homosapiens_catalogue_mh.py) is used to automatically retrieve the latest HLA data from the github repository above, and the data retrieved is transformed into json form, where the data is encoded as:
    ```
    {
        gene name: {
            allele num 1: {
                allele num 2: {
                    possibly allele num 3: {
                        possibly allele num 4: {}
                    }
                }
            }
        }
    }
    ```

## Mus musculus

### Valid TR genes

- Data source: IMGT/GENE-DB (https://www.imgt.org/genedb/)
    - Use interactive GENE-DB interface to retrieve all Mus musculus TR genes
    - Then, in order to get a list of all known alleles of each gene, request the nucleotide sequence of all functional, ORF and pseudogene alleles in FASTA format
    - Downloaded FASTA data found at: https://github.com/yutanagano/tidytcells/tree/main/data/musmusculus_tr.fasta
- Last date of access: 22nd September 2023
- Data parsing:
    - Using the FASTA data retrieved as above, use the data processing script (found at: https://github.com/yutanagano/tidytcells/tree/main/scripts/musmusculus_catalogue_tr.py) to transform the data to json form, where the data is encoded as:
    ```
    {
        gene name: {
            allele number: functionality
        }
    }
    ```

### Valid MH genes

- Data source: IMGT Repertoire (https://www.imgt.org/IMGTrepertoireMH/index.php?section=LocusGenes&repertoire=Chromosomal%20localization&species=mouse)
    - Tabular data from the above website is copied onto an ODS spreadsheet (found at: https://github.com/yutanagano/tidytcells/tree/main/data/musmusculus_mh.ods)
- Last date of access: 22nd September 2023
- Data parsing:
    - The data processing script (found at: https://github.com/yutanagano/tidytcells/tree/main/scripts/musmusculus_catalogue_mh.py) transforms the data from the spreadsheet as above into a json, where the data is encoded as:
    ```
    [
        gene name
    ]
    ```