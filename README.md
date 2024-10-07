<div align="center">

<img src="https://raw.githubusercontent.com/yutanagano/tidytcells/main/tidytcells.png" width=700>
<br><br>

![Tests](https://github.com/yutanagano/tidytcells/actions/workflows/tests.yaml/badge.svg)
[![Docs](https://readthedocs.org/projects/tidytcells/badge/?version=latest)](https://tidytcells.readthedocs.io)
[![License](https://img.shields.io/badge/license-MIT-blue)](https://github.com/yutanagano/tidytcells?tab=MIT-1-ov-file#readme)
[![DOI](https://img.shields.io/badge/DOI-10.3389/fimmu.2023.1276106-pink)](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2023.1276106)

### Check out the [documentation page](https://tidytcells.readthedocs.io).

</div>

---

`tidytcells` is a lightweight python package that cleans and standardizes T cell receptor (TR) and Major Histocompatibility (MH) data to be [IMGT](https://www.imgt.org/)-compliant.
The main purpose of the package is to solve the problem of parsing and collating together non-standardized TR datasets.
It is often difficult to compile TR data from multiple sources because the formats/nomenclature of how each dataset encodes TR and MH gene names are slightly different, or even inconsistent within themselves.
`tidytcells` can ameliorate this issue by auto-correcting and auto-standardizing your data.

## Installation

```bash
$ pip install tidytcells
```

## Citing tidytcells

Please cite [our manuscript](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2023.1276106).

### BibTex
```bibtex
@ARTICLE{10.3389/fimmu.2023.1276106,
         AUTHOR={Nagano, Yuta  and Chain, Benjamin },
         TITLE={tidytcells: standardizer for TR/MH nomenclature},
         JOURNAL={Frontiers in Immunology},
         VOLUME={14},
         YEAR={2023},
         URL={https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2023.1276106},
         DOI={10.3389/fimmu.2023.1276106},
         ISSN={1664-3224}
}
```
