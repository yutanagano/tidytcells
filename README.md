<h1 align="center">
    <img src="https://raw.githubusercontent.com/yutanagano/tidytcells/main/tidytcells.png" width=700>
</h1>

![Tests](https://github.com/yutanagano/tidytcells/actions/workflows/tests.yaml/badge.svg)
[![Docs](https://readthedocs.org/projects/tidytcells/badge/?version=latest)](https://tidytcells.readthedocs.io)
[![License](https://img.shields.io/badge/license-MIT-blue)](https://github.com/yutanagano/tidytcells?tab=MIT-1-ov-file#readme)

`tidytcells` is a lightweight python package that cleans and standardizes T cell receptor (TR) and Major Histocompatibility (MH) data to be [IMGT](https://www.imgt.org/)-compliant.
The main purpose of the package is to solve the problem of parsing and collating together non-standardized TR datasets.
It is often difficult to compile TR data from multiple sources because the formats/nomenclature of how each dataset encodes TR and MH gene names are slightly different, or even inconsistent within themselves.
`tidytcells` can ameliorate this issue by auto-correcting and auto-standardizing your data!
Check out the [documentation page](https://tidytcells.readthedocs.io).

## Installation

### Via [PyPI](https://pypi.org/project/tidytcells/) (recommended)

`tidytcells` can be installed using `pip`:

```bash
$ pip install tidytcells
```

### From [source](https://github.com/yutanagano/tidytcells)

The source code for the package is available [on Github](https://github.com/yutanagano/tidytcells).
To install from source, clone the git repository, and run:

```bash
$ pip install .
```

from inside the project root directory.

## Useful links

- [Documentation](https://tidytcells.readthedocs.io)
- [PyPI page](https://pypi.org/project/tidytcells)
