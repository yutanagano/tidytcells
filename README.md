# tidytcells

> DISCLAIMER: This package currently only supports parsing of human TCR and MHC
> gene data. Support for more species is planned for the future.

`tidytcells` is a lightweight Python package written for bioinformaticians who
work with T cell receptor (TCR) data. The main purpose of the package is to
solve the problem of parsing and collating together non-standardised TCR
datasets. It is often difficult to compile TCR data from multiple sources
because the formats/nomenclature of how each dataset encodes TCR and MHC gene
names are slightly different, or even inconsistent within themselves.
`tidytcells` attempts to ameliorate this issue by providing simple functions
that can, among other things, standardise TCR and MHC gene names to be
IMGT-compliant.

The package is currently in its alpha stage in development. More thorough
documentation will follow soon.

## Table of contents

1. [Installation](#installation)
1. [API reference](#api-reference)
    - [`mhc` submodule](#mhc-submodule)
    - [`tcr` submodule](#tcr-submodule)
2. [Example usage](#example-usage)
3. [Contributions](#contributions)

## Installation

|Installation method|Comments|
|-|-|
|[PyPI](https://pypi.org/project/tidytcells/) (recommended)|`pip install tidytcells`|
|[source](https://github.com/yutanagano/tidytcells)|Download the source code and run `pip install .` from the project root directory|

## API reference

`tidytcells` currently comes with two submodules that can directly be accessed
from the parent module: `tidytcells.mhc` and `tidytcells.tcr`.

### `.mhc` submodule

#### `.standardise(gene_name: str, species: str)`

Attempt to standardise `gene_name` for `species` to be IMGT-compliant.

|Parameters|Type|Value|
|-|-|-|
|`gene_name`|`str`|Potentially non-standard name for MHC gene|
|`species`|`str`|Species to which the MHC gene belongs (see [below](#species-names))|

|Return type|Value|
|-|-|
|`Tuple[Union[str, None], Union[str, None]]`|If the specified `species` is supported, and `gene_name` could be standardised, then return a tuple containing the standardised gene name decomposed into two parts: 1) the name of the gene specific to the level of the protein, and 2) (if any) further valid specifier fields. If `species` is unsupported, then return a tuple with the `gene_name` as is for the first element, and `None` for the second element. Else return the tuple `(None, None)`. See [example usage](#example-usage).|

#### `.get_chain(gene_name: str)`

> NOTE: This function currently only supports HLA gene names.

Given an IMGT-compliant MHC gene name `gene_name`, detect whether it codes for
an alpha chain or a beta chain.

|Parameters|Type|Value|
|-|-|-|
|`gene_name`|`str`|IMGT-compliant MHC gene name|

|Return type|Value|
|-|-|
|`Union[str, None]`|`'alpha'` or `'beta'` if `gene_name` recognised, else `None`|

#### `.classify(gene_name: str)`

> NOTE: This function currently only supports HLA gene names.

Given an IMGT-compliant MHC gene name `gene_name`, detect whether it comprises
a class I or II MHC receptor.

|Parameters|Type|Value|
|-|-|-|
|`gene_name`|`str`|IMGT-compliant MHC gene name|

|Return type|Value|
|-|-|
|`Union[int, None]`|`1` or `2` if `gene_name` recognised, else `None`|

### `.tcr` submodule

#### `.standardise(gene_name: str, species: str)`

Attempt to standardise `gene_name` for `species` to be IMGT-compliant.

|Parameters|Type|Value|
|-|-|-|
|`gene_name`|`str`|Potentially non-standard name for TCR gene|
|`species`|`str`|Species to which the TCR gene belongs (see [below](#species-names))|

|Return type|Value|
|-|-|
|`Union[str, None]`|If the specified `species` is supported, and `gene_name` could be standardised, then return the standardised gene name. If `species` is unsupported, then return `gene_name` as is. Else return `None`.|

### Species names

For all functions that expect a species to be specified via a string, the
species should be referred to by its binomial name (genus followed by species),
CamelCased, with no space between the two parts (e.g. `'HomoSapiens'`).

## Example usage

```Python
import tidytcells


# --- MHC parsing ---


tidytcells.mhc.standardise('HLA-A', 'HomoSapiens')
# > ('HLA-A', None)

tidytcells.mhc.standardise('B07', 'HomoSapiens')
# > ('HLA-B*07', None)

tidytcells.mhc.standardise('DRA*01:01:01', 'HomoSapiens')
# > ('HLA-DRA*01:01', ':01')

tidycells.mhc.get_chain('HLA-A')
# > 'alpha'

tidycells.mhc.classify('HLA-DRB1*01:01')
# > 2


# --- TCR parsing ---


tidycells.tcr.standardise('TCRAV32S1', 'HomoSapiens')
# > 'TRAV25'
```

## Contributions

Please feel free to contribute by submitting bug reports and pull requests.