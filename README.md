# tidytcells

> NOTE: This package is currently in the alpha stage of development.

> NOTE: This package currently only supports parsing of human TCR and MHC gene
> data. Support for more species is planned for the future.

`tidytcells` is a lightweight Python package written for bioinformaticians who
work with T cell receptor (TCR) data. The main purpose of the package is to
solve the problem of parsing and collating together non-standardised TCR
datasets. It is often difficult to compile TCR data from multiple sources
because the formats/nomenclature of how each dataset encodes TCR and MHC gene
names are slightly different, or even inconsistent within themselves.
`tidytcells` attempts to ameliorate this issue by providing simple functions
that can, among other things, standardise TCR and MHC gene names to be
IMGT-compliant.

## Useful links

- [PyPI page](https://pypi.org/project/tidytcells/)