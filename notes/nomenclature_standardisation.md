# Nomenclature standardisation: implementation notes

## TCR genes

### Human TCR genes

#### Dealing with deprecated gene names

In order to resolve TCR gene names specified via deprecated nomenclature (which include ambiguous deprecated symbols which could refer to one or more genes), we use data on currently accepted TCR gene symbols and currently/previously accepted symbols from [HGNC](www.genenames.org).

We downloaded HGNC data via their [custom downloads](https://www.genenames.org/download/custom/) web service, with the following column data:

* HGNC ID
* Approved symbol
* Approved name
* Status
* Locus type
* Previous symbols
* Alias symbols
* Gene group name

where data on withdrawn genes were excluded from the download.

This allows you to download the specified column data for every entry in HGNC as a TSV (tab-separated values) file.
This source table was parsed (using code in `notebooks/homosapiens_catalogue_tcr.ipynb`) to create a dictionary of deprecated/alias TCR gene symbols and their modern translation.

#### Cross-checking standardisation attempt with list of known TCR genes

As a final step of TCR gene nomenclature standardisation, we make sure that the attempted standardisation of a TCR gene name actually maps to an existing gene.
In order to ensure this fact, the standardisation programme cross-checks its attempted fix of a TCR gene with a table of known, existing TCR genes (+ their alleles).
If an attempted fix does not look like any real TCR gene/allele, it considers the input to be nonsensical and will throw it away.

The data listing all known TCR genes and alleles is sourced from [IMGT](www.imgt.org).
The raw source data is parsed by `notebooks/homosapiens_catalogue_tcr` to create reference jsons for the standardiser programme.

### Mus musculus TCR genes

Same as humans, but we could not find a good source for identifying deprecated names or aliases.

## MHC Genes

### Human MHC (HLA) genes

#### Dealing with deprecated gene names

Similarly to the case with human TCR genes, we used downloaded data from [HGNC](www.genenames.org) to resolve deprecated gene names to the current one.

#### Cross-checking standardisation attempt with list of known HLA genes/groups

Similarly to the case with TCR genes, we needed a set of reference data of approved symbols for known HLA genes.
This data was downloaded from the [IPD-IMGT](https://www.ebi.ac.uk/ipd/imgt/hla/download/) database service and parsed by `notebooks/homosapiens_catalogue_mhc` to create a hierarchical tree of possible HLA names.