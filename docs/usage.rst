Usage
=====

Standardising TCR/MHC gene names
--------------------------------

These are :py:mod:`tidytcells`' primary usecases.
The API for standardising TCR and MHC genes are similar: :py:func:`tidytcells.tcr.standardise` and :py:func:`tidytcells.mhc.standardise` respectively.
At their core, these functions take as input a potentially non-standard string representation of a gene or allele and returns (if possible) an IMGT-standardised version of the gene/allele name.
In both cases, '``standardize``' is a valid alias of the function as well.

Other MHC utilities
-------------------

:py:mod:`tidytcells` also provides some extra utilities in the :py:mod:`mhc <tidytcells.mhc>` module, including :py:func:`get_chain <tidytcells.mhc.get_chain>` and :py:func:`get_class <tidytcells.mhc.get_class>`.

.. _example_usage:

Example usage
-------------

>>> import tidytcells
>>> # --- TCR parsing ---
>>> tidytcells.tcr.standardise('TCRAV32S1', 'homosapiens')
'TRAV25'
>>> tidytcells.tcr.standardise('TRBV1*01', 'homosapiens', enforce_functional=True)
None
>>> tidytcells.tcr.standardise('TRAJ12*01', 'musmusculus', precision='gene')
'TRAJ12'
>>> # --- MHC parsing ---
>>> tidytcells.mhc.standardise('HLA-A', 'homosapiens')
'HLA-A'
>>> tidytcells.mhc.standardise('A1', 'homosapiens')
'HLA-A*01'
>>> tidytcells.mhc.standardise('HLA-B*07:02:01:01', 'homosapiens', precision='protein')
('HLA-B*07:02', ':01:01')
>>> tidytcells.mhc.standardise('HLA-DR1BL')
'HLA-DRB9'
>>> tidytcells.mhc.get_chain('HLA-A')
'alpha'
>>> tidytcells.mhc.get_class('HLA-DRB1*01:01')
2

.. _supported_species:

Supported species and species strings
-------------------------------------

For all functions that expect a species to be specified via a string, the species should be referred to by its `binomial name <https://en.wikipedia.org/wiki/Binomial_nomenclature>`_, lowercased, with no space between the two parts (e.g. ``'homosapiens'``).
That being said, for backwards-compatibility and user friendliness, any combination of camelcasing, uppercasing and the use of whitespace is technically allowed and recognised by the software.

Below is a list of currently supported species:

:py:func:`tidytcells.mhc.get_chain`
    - ``homosapiens``

:py:func:`tidytcells.mhc.get_class`
    - ``homosapiens``

:py:func:`tidytcells.mhc.standardise`
    - ``homosapiens``
    - ``musmusculus``

:py:func:`tidytcells.tcr.standardise`
    - ``homosapiens``
    - ``musmusculus``