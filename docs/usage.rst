Usage
=====

``tidytcells``' structure
-------------------------

:py:mod:`tidytcells` is comprised of the :py:mod:`tcr <tidytcells.tcr>` and :py:mod:`mhc <tidytcells.mhc>` submodules.
Each module provides useful functions that can help you standardize or otherwise deal with data from their respective category.
For ease of use, function APIs are standardised accross modules wherever possible- for example, each module has a function named ``standardise`` (see below) which standardises data from each category to be `IMGT <https://www.imgt.org/>`_-compliant.
Refer to :ref:`here <api>` for a full review of :py:mod:`tidytcells`' API.

Standardising TCR/MHC gene names
--------------------------------

This is :py:mod:`tidytcells`' primary usecase.
Both the :py:func:`tcr <tidytcells.tcr.standardise>` and :py:func:`mhc <tidytcells.mhc.standardise>` module provides a function named ``standardise`` which take as input a potentially non-standard string representation of a gene or allele and returns (if possible) an `IMGT <https://www.imgt.org/>`_-standardised version of the gene/allele name.
In both cases, '``standardize``' is a valid alias of the function as well.

Querying from `IMGT <https://www.imgt.org/>`_ TCR/MHC genes or alleles
----------------------------------------------------------------------

:py:mod:`tidytcells` also provides the nifty functions :py:func:`tidytcells.tcr.query` and :py:func:`tidytcells.mhc.query` that allows users to obtain a list (actually a ``FrozenSet``) of `IMGT <https://www.imgt.org/>`_ gene/allele names from the respective categories.
The functions allow the user to provide various constraints relating to the genes/alleles' functionalities and names to filter the query results as well.
The ``query`` functions can be useful when checking if a particular dataset covers all the TCR or MHC genes, or counting how many genes fulfill a particular set of constraints.

Other MHC utilities
-------------------

The :py:mod:`mhc <tidytcells.mhc>` module provides a couple more extra goodies, including :py:func:`get_chain <tidytcells.mhc.get_chain>` and :py:func:`get_class <tidytcells.mhc.get_class>`, each with self-explanatory names.

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