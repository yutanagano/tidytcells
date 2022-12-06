# Nomenclature standardisation: implementation notes

## TCR genes

### Human TCR genes

#### Dealing with deprecated gene names

In order to resolve TCR gene names specified via deprecated nomenclature (which
include ambiguous deprecated symbols which could refer to one or more genes),
we use data on currently accepted TCR gene symbols and currently/previously
accepted symbols from [HGNC](www.genenames.org).

We downloaded HGNC data via their
[custom downloads](https://www.genenames.org/download/custom/) web service,
with the following column data:

* HGNC ID
* Approved symbol
* Approved name
* Status
* Locus type
* Previous symbols
* Alias symbols
* Gene group name

where data on withdrawn genes were excluded from the download.

This allows you to download the specified column data for every entry in HGNC
as a TSV (tab-separated values) file. This source table was parsed (using code
in `notebooks/catalogue_human_tcr_gene_data.ipynb`) to create a dictionary of
deprecated/alias TCR gene symbols and their modern translation. Cases of
ambiguous naming were handled as documented in the main README.

#### Cross-checking standardisation attempt with list of known TCR genes

As a final step of TCR gene nomenclature standardisation, we make sure that the
attempted standardisation of a TCR gene name actually maps to an existing gene.
In order to ensure this fact, the standardisation programme cross-checks its
attempted fix of a TCR gene with a table of known, existing TCR genes (+ their
alleles). If an attempted fix does not look like any real TCR gene/allele, it
considers the input to be nonsensical and will throw it away.

The data listing all known TCR genes and alleles is sourced from
[IMGT](www.imgt.org). The raw source data is parsed by
`notebooks/catalogue_human_tcr_gene_data` to create reference tables for the
standardiser programme.

## MHC Genes

### Human MHC (HLA) genes

#### Dealing with deprecated gene names

Similarly to the case with human TCR genes, we used downloaded data from
[HGNC](www.genenames.org) to resolve deprecated gene names to the current one.

#### Cross-checking standardisation attempt with list of known HLA genes/groups

Similarly to the case with TCR genes, we needed a set of reference data of
approved symbols for known HLA genes. This data was downloaded from the
[IPD-IMGT](https://www.ebi.ac.uk/ipd/imgt/hla/download/) database service and
parsed by `notebooks/catalogue_human_mhc_gene_data` to create a hierarchical
tree of possible HLA names.

Furthermore, in order to get a list of all known G and P groups of HLA alleles,
we downloaded an XML with G and P group information and parsed it using
`notebooks/catalogue_human_mhc_gene_data` to generate a list of all G and P
groups.