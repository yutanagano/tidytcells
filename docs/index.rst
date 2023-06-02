.. tidytcells documentation master file, created by
   sphinx-quickstart on Sun Nov  6 16:32:10 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to tidytcells's documentation!
======================================

:py:mod:`tidytcells` is a lightweight python package that cleans and standardizes T cell receptor (TCR) and Major Histocompatibility Complex (MHC) data to be `IMGT <https://www.imgt.org/>`_-compliant.
The main purpose of the package is to solve the problem of parsing and collating together non-standardized TCR datasets.
It is often difficult to compile TCR data from multiple sources because the formats/nomenclature of how each dataset encodes TCR and MHC gene names are slightly different, or even inconsistent within themselves.
:py:mod:`tidytcells` can ameliorate this issue by auto-correcting and auto-standardizing your data!

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