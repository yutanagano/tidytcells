.. tidytcells documentation master file, created by
   sphinx-quickstart on Sun Nov  6 16:32:10 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to tidytcells's documentation!
======================================

.. note::

   This package is currently in the alpha stage of development.

.. note::

   Support for species other than Homo sapiens is currently limited
   (see :ref:`_supported_species`). Support for more species is planned for the
   future.

``tidytcells`` is a lightweight Python package written for bioinformaticians
who work with T cell receptor (TCR) data. The main purpose of the package is to
solve the problem of parsing and collating together non-standardised TCR
datasets. It is often difficult to compile TCR data from multiple sources
because the formats/nomenclature of how each dataset encodes TCR and MHC gene
names are slightly different, or even inconsistent within themselves.
``tidytcells`` attempts to ameliorate this issue by providing simple functions
that can standardise TCR and MHC gene symbols to their officially accepted
versions, as defined by `IMGT <https://www.imgt.org/>`_,
`HGNC <https://www.genenames.org/>`_, or other authorities on gene
nomenclature.

Contents
--------

.. toctree::
   :maxdepth: 2
   
   installation
   usage
   api
   contribute

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`