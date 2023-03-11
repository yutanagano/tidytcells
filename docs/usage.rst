Usage
=====

:py:mod:`tidytcells`' structure
-------------------------

:py:mod:`tidytcells` is comprised of several modules, each of which provide a set of functions that help process a particular type of data that bioinformaticians working on TCR data may come accross.

The submodules are:

+-------------------------------+----------------------------------------------------------+
| Submodule                     | For                                                      |
+===============================+==========================================================+
| :py:mod:`tidytcells.aa`       | General amino acid sequence data (e.g. peptide epitopes) |
+-------------------------------+----------------------------------------------------------+
| :py:mod:`tidytcells.junction` | TCR junction (CDR3) amino acid data                      |
+-------------------------------+----------------------------------------------------------+
| :py:mod:`tidytcells.mhc`      | MHC gene/allele data                                     |
+-------------------------------+----------------------------------------------------------+
| :py:mod:`tidytcells.tcr`      | TCR gene/allele data                                     |
+-------------------------------+----------------------------------------------------------+

For ease of use, function APIs are standardised accross modules wherever possible- for example, each module has a function named ``standardise`` (see below) which standardises data from each category to be `IMGT <https://www.imgt.org/>`_-compliant.
Refer to :ref:`here <api>` for a full review of :py:mod:`tidytcells`' API.

Standardising TCR/MHC data using :py:mod:`tidytcells` and `pandas <https://pandas.pydata.org/>`_
------------------------------------------------------------------------------------------------

This is :py:mod:`tidytcells`' primary usecase.

Since each of :py:mod:`tidytcells`' submodules provide a ``standardise`` (``standardize`` is a valid alias as well) function that automates data cleaning in their respective data category, these functions can be used in ensemble to clean a whole dataset of TCR/MHC data.
Now, these ``standardise`` functions can be used on their own to clean individual pieces of data- that is for example:

>>> import tidytcells as tt
>>> orig = "A1"
>>> cleaned = tt.mhc.standardise(orig)
>>> cleaned
'HLA-A*01'

However, in real-life scenarios one would like to clean a whole set of data contained in a table.
This can be achieved in a fairly straightforward manner by using :py:mod:`tidytcells` in conjunction with a data analysis tool like `pandas <https://pandas.pydata.org/>`_.
Pandas provides a nice way to blanket-apply data transformation functions to multiple ``DataFrame`` cells through their ``Series.map`` and ``DataFrame.applymap`` methods.
Therefore, given a table of TCR and MHC data:

>>> import pandas as pd
>>> df = pd.DataFrame(
...     data=[
...         ["TRBV13*01",    "CASSYLPGQGDHYSNQPQHF", "trbj1-5*01"],
...         ["TCRBV28S1*01", "CASSLGQSGANVLTF",      "TRBJ2-6*01"],
...         ["unknown",      "ASSDWGSQNTLY",         "TRBJ2-4*01"]
...     ],
...     columns=["v", "junction", "j"]
... )
>>> df
              v              junction           j
0     TRBV13*01  CASSYLPGQGDHYSNQPQHF  trbj1-5*01
1  TCRBV28S1*01       CASSLGQSGANVLTF  TRBJ2-6*01
2       unknown          ASSDWGSQNTLY  TRBJ2-4*01

One can apply the ``standardise`` functions from :py:mod:`tidytcells` over the whole table at once, like so:

>>> cleaned = df.copy()
>>> cleaned[["v", "j"]] = df[["v", "j"]].applymap(tt.tcr.standardise)
>>> cleaned["junction"] = df["junction"].map(tt.junction.standardise)
>>> cleaned
           v              junction           j
0  TRBV13*01  CASSYLPGQGDHYSNQPQHF  TRBJ1-5*01
1  TRBV28*01       CASSLGQSGANVLTF  TRBJ2-6*01
2       None        CASSDWGSQNTLYF  TRBJ2-4*01

To apply the functions with optional arguments, one can wrap the ``standardise`` functions using lambda functions (see below).
For use cases that require more flexibility, one could even define a wrapper function explicitly in the code.

>>> cleaned = df.copy()
>>> cleaned[["v", "j"]] = df[["v", "j"]].applymap(
...     lambda x: tt.tcr.standardise(
...         gene=x,
...         species="homosapiens",
...         precision="gene"
...     )
... )
>>> cleaned["junction"] = df["junction"].map(
...     lambda x: tt.junction.standardise(
...         seq=x,
...         strict=True
...     )
... )
>>> cleaned
        v              junction        j
0  TRBV13  CASSYLPGQGDHYSNQPQHF  TRBJ1-5
1  TRBV28       CASSLGQSGANVLTF  TRBJ2-6
2    None                  None  TRBJ2-4

For more complete documentations of the ``standardise`` functions, refer to :ref:`the api reference <api>`.

Querying from `IMGT <https://www.imgt.org/>`_ TCR/MHC genes or alleles
----------------------------------------------------------------------

:py:mod:`tidytcells` also provides the nifty functions :py:func:`tidytcells.tcr.query` and :py:func:`tidytcells.mhc.query` that allows users to obtain a list (actually a ``FrozenSet``) of `IMGT <https://www.imgt.org/>`_ gene/allele names from the respective categories.
The functions allow the user to provide various constraints relating to the genes/alleles' functionalities and names to filter the query results as well.
The ``query`` functions can be useful when checking if a particular dataset covers all the TCR or MHC genes, or counting how many genes fulfill a particular set of constraints.

Other MHC utilities
-------------------

The :py:mod:`mhc <tidytcells.mhc>` module provides a couple more extra goodies, including :py:func:`get_chain <tidytcells.mhc.get_chain>` and :py:func:`get_class <tidytcells.mhc.get_class>`, each with self-explanatory names.