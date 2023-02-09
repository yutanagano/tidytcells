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
>>> tidytcells.tcr.standardise('TCRAV32S1', 'HomoSapiens')
'TRAV25'
>>> # --- MHC parsing ---
>>> tidytcells.mhc.standardise('HLA-A', 'HomoSapiens')
'HLA-A'
>>> tidytcells.mhc.standardise('HLA-DR1BL', 'HomoSapiens')
'HLA-DRB9'
>>> tidycells.mhc.get_chain('HLA-A')
'alpha'
>>> tidycells.mhc.get_class('HLA-DRB1*01:01')
2

.. _supported_species:

Supported species and species strings
-------------------------------------

For all functions that expect a species to be specified via a string, the species should be referred to by its `binomial name <https://en.wikipedia.org/wiki/Binomial_nomenclature>`_, `CamelCased <https://en.wikipedia.org/wiki/Camel_case>`_ (with the first character capitalised), with no space between the two parts (e.g. ``'HomoSapiens'``).

Below is a list of currently supported species:

:py:func:`tidytcells.mhc.classify`
    - ``HomoSapiens``

:py:func:`tidytcells.mhc.get_chain`
    - ``HomoSapiens``

:py:func:`tidytcells.mhc.standardise`
    - ``HomoSapiens``

:py:func:`tidytcells.tcr.standardise`
    - ``HomoSapiens``
    - ``MusMusculus``