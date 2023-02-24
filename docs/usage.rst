Usage
=====

.. note::
    
    As stated on the home page, there is currently limited support for non-human species.
    Support for more species are planned for the future.

Standardising TCR/MHC gene names
--------------------------------

``tidytcells`` provides a similar API for standardising TCR and MHC genes.
The function for standardising TCR and MHC genes are ``tidytcells.tcr.standardise`` and ``tidytcells.mhc.standardise`` respectively.

.. autofunction:: tidytcells.tcr.standardise

.. autofunction:: tidytcells.mhc.standardise

Other MHC utilities
-------------------

.. autofunction:: tidytcells.mhc.get_chain

.. autofunction:: tidytcells.mhc.get_class

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

Below is a list of currently supported species:

:py:func:`tidytcells.mhc.classify`
    - ``homosapiens``

:py:func:`tidytcells.mhc.get_chain`
    - ``homosapiens``

:py:func:`tidytcells.mhc.standardise`
    - ``homosapiens``

:py:func:`tidytcells.tcr.standardise`
    - ``homosapiens``
    - ``musmusculus``