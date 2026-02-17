Usage
=====

:py:mod:`tidytcells`' structure
-------------------------------

:py:mod:`tidytcells` is comprised of several modules, each of which provide a set of functions that help process a particular type of data that bioinformaticians working on T cell receptor (TR) or Major Histocompatibility (MH) data may come accross.

The submodules are:

+-------------------------------+----------------------------------------------------------+
| Submodule                     | For                                                      |
+===============================+==========================================================+
| :py:mod:`tidytcells.tr`       | TR gene/allele data                                      |
+-------------------------------+----------------------------------------------------------+
| :py:mod:`tidytcells.ig`       | IG gene/allele data                                      |
+-------------------------------+----------------------------------------------------------+
| :py:mod:`tidytcells.junction` | TR JUNCTION or CDR3-IMGT amino acid sequence data        |
+-------------------------------+----------------------------------------------------------+
| :py:mod:`tidytcells.mh`       | MH gene/allele data                                      |
+-------------------------------+----------------------------------------------------------+

.. tip::
   
   The :py:mod:`tidytcells.ig` submodule is newly added! It provides functionality for standardizing, querying, and retrieving amino acid sequences for immunoglobulin genes/alleles, similar to the existing TR and MH modules. Thanks to `Lonneke <https://github.com/LonnekeScheffer>`_ for implementing this module!

For ease of use, function APIs are standardized accross modules wherever possible- for example, each module has a function named ``standardize`` (see below) which standardizes data from each category to be IMGT-compliant (`IMGT/GENE-DB <https://www.imgt.org/genedb/>`_, `IMGT Repertoire <https://www.imgt.org/IMGTrepertoire/>`_).
Refer to :ref:`here <api>` for a full review of :py:mod:`tidytcells`' API.

Standardizing TR/junction/peptide-MH data using :py:mod:`tidytcells` and `pandas <https://pandas.pydata.org/>`_
---------------------------------------------------------------------------------------------------------------

This is :py:mod:`tidytcells`' primary usecase.

Note that :py:mod:`tidytcells` also provides functions to standardize IG genes/alleles in the same way as described below. However, for the sake of simplicity, this example focuses on TR data.

Since each of :py:mod:`tidytcells`' submodules provide a ``standardize`` (``standardise`` is a valid alias as well) function that automates data cleaning in their respective data category, these functions can be used in ensemble to clean a whole dataset of TR/MH data.
Now, these ``standardize`` functions can be used on their own to clean individual pieces of data- that is for example:

>>> import tidytcells as tt
>>> orig = "A1"
>>> cleaned = tt.mh.standardize(orig)
>>> cleaned.symbol
'HLA-A*01'

The result of each ``standardize`` function is a wrapper object which will have various properties and utility
methods for retrieving further information. Please see the :py:class:`~tidytcells.result.ReceptorGene`,
:py:class:`~tidytcells.result.Junction` and :py:class:`~tidytcells.result.MhGene` classes for more detailed documentation.



In real-life scenarios one would like to clean a whole set of data contained in a table.
This can be achieved in a fairly straightforward manner by using :py:mod:`tidytcells` in conjunction with a data analysis tool like `pandas <https://pandas.pydata.org/>`_.
Pandas provides a nice way to blanket-apply data transformation functions to multiple ``DataFrame`` cells through their ``Series.map`` and ``DataFrame.map`` methods.
For example, given a table of TR/junction data (a similar procedure would work for tables with peptide-MH data as well):

>>> import pandas as pd
>>> df = pd.DataFrame(
...     data=[
...         ["TRBV13",    "CASSYLPGQGDHYSNQPQHF", "trbj1-5*01"],
...         ["TCRBV28S1*01", "CASSLGQSGANVLTF",      "TRBJ2-6*01"],
...         ["unknown",      "ASSDWGSQNTLY",         "TRBJ2-4*01"]
...     ],
...     columns=["v_orig", "junction_orig", "j_orig"]
... )
>>> df
         v_orig         junction_orig      j_orig
0        TRBV13  CASSYLPGQGDHYSNQPQHF  trbj1-5*01
1  TCRBV28S1*01       CASSLGQSGANVLTF  TRBJ2-6*01
2       unknown          ASSDWGSQNTLY  TRBJ2-4*01

One can apply the ``standardize`` functions from :py:mod:`tidytcells` over the whole table at once.
Use 'map' for standardization with default parameters, and 'apply' when parameters need to be set:

>>> cleaned = df.copy()
>>> cleaned[["v", "j"]] = df[["v_orig", "j_orig"]].map(tt.tr.standardize)
>>> cleaned["junction"] = df["cdr3_orig"].apply(junction.standardize, locus="TRB")
>>> cleaned[["v", "junction", "j"]]
           v              junction           j
0     TRBV13  CASSYLPGQGDHYSNQPQHF  TRBJ1-5*01
1  TRBV28*01       CASSLGQSGANVLTF  TRBJ2-6*01
2                   CASSDWGSQNTLYF  TRBJ2-4*01

Alternatively, for full control over the input parameters and retrieving different output values from the results objects,
one one can wrap the ``standardize`` functions using lambda functions (see below).
For use cases that require more flexibility, one could even define a wrapper function explicitly in the code.

>>> cleaned = df.copy()
>>> cleaned["v_allele"] = cleaned["v_orig"].map(
...     lambda x: tt.tr.standardize(
...         symbol=x,
...     ).allele
... )
>>> cleaned["v_gene"] = cleaned["v_orig"].map(
...     lambda x: tt.tr.standardize(
...         symbol=x,
...     ).gene
... )
>>> cleaned["v_subgroup"] = cleaned["v_orig"].map(
...     lambda x: tt.tr.standardize(
...         symbol=x,
...     ).subgroup
... )
>>> cleaned[["v_orig", "v_allele", "v_gene", "v_subgroup"]]
         v_orig   v_allele  v_gene v_subgroup
0        TRBV13       None  TRBV13     TRBV13
1  TCRBV28S1*01  TRBV28*01  TRBV28     TRBV28
2       unknown       None    None       None

>>> cleaned["cdr3"] = cleaned["junction_orig"].map(
...     lambda x: tt.junction.standardize(
...         seq=x,
...         locus="TRB",
...     ).cdr3
... )
>>> cleaned["junction"] = cleaned["junction_orig"].map(
...     lambda x: tt.junction.standardize(
...         seq=x,
...         locus="TRB",
...     ).cdr3
... )
>>> cleaned["junction_orig", "junction", "cdr3"]
          junction_orig              junction                cdr3
0  CASSYLPGQGDHYSNQPQHF  CASSYLPGQGDHYSNQPQHF  ASSYLPGQGDHYSNQPQH
1       CASSLGQSGANVLTF       CASSLGQSGANVLTF       ASSLGQSGANVLT
2          ASSDWGSQNTLY        CASSDWGSQNTLYF        ASSDWGSQNTLY


For more complete documentations of the ``standardize`` functions, refer to :ref:`the api reference <api>`.

Querying from `IMGT TR/MH/IG genes or alleles <https://www.imgt.org/IMGTrepertoire/>`_
--------------------------------------------------------------------------------------

:py:mod:`tidytcells` also provides the nifty functions :py:func:`tidytcells.tr.query`, :py:func:`tidytcells.mh.query`, and :py:func:`tidytcells.ig.query` that allows users to obtain a list (actually a ``FrozenSet``) of `IMGT gene/allele names <https://www.imgt.org/IMGTrepertoire/>`_ from the respective categories.
The functions allow the user to provide various constraints relating to the genes/alleles' functionalities and names to filter the query results as well.
The ``query`` functions can be useful when checking if a particular dataset covers all the TR, MH, or IG genes, or counting how many genes fulfill a particular set of constraints.
Since :py:mod:`tidytcells` has a local copy of all relevant data pulled directly from `IMGT's GENE-DB <https://www.imgt.org/genedb/>`_ (and updated with every new release), queries are blazingly fast and do not require an internet connection.




..
    To do: rewrite these parts to match with new results objects / MRO

    Querying TR/IG gene amino acid sequence data from `IMGT GENE-DB <https://www.imgt.org/genedb/>`_
    -------------------------------------------------------------------------------------------------

    Sometimes, you have a T cell receptor or immunoglobulin represented as its V and J gene usages and its junction sequences, but you want to represent it in terms of its amino acid sequence.
    In such situations, the :py:func:`tidytcells.tr.get_aa_sequence` and :py:func:`tidytcells.ig.get_aa_sequence` functions can help.
    These functions allow you to query amino acid sequence data for any functional TR or IG gene.
    The functions provide sequence data for the whole gene exome, as well as certain important regions (e.g. CDR1 and CDR2 in the V genes).
    The data is pulled from IMGT's `GENE-DB <https://www.imgt.org/genedb/>`_, and as is with the case with the :py:func:`tidytcells.tr.query`, :py:func:`tidytcells.mh.query`, and :py:func:`tidytcells.ig.query`, all relevant data exists locally within :py:mod:`tidytcells` (and updated with every new release), so the queries are blazingly fast and requires no internet connection.

    Other MH utilities
    ------------------

    The :py:mod:`mh <tidytcells.mh>` module provides a couple more extra goodies, including :py:func:`get_chain <tidytcells.mh.get_chain>` and :py:func:`get_class <tidytcells.mh.get_class>`, each with self-explanatory names.