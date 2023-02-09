tidytcells
==========

.. note::

    This package is currently in the alpha stage of development.

.. note::

    Support for species other than Homo sapiens is currently limited (see `the docs <https://tidytcells.readthedocs.io>`_).
    Support for more species is planned for the future.

``tidytcells`` is a lightweight Python package written for bioinformaticians who work with T cell receptor (TCR) data.
The main purpose of the package is to solve the problem of parsing and collating together non-standardised TCR datasets.
It is often difficult to compile TCR data from multiple sources because the formats/nomenclature of how each dataset encodes TCR and MHC gene names are slightly different, or even inconsistent within themselves.
``tidytcells`` attempts to ameliorate this issue by providing simple functions that can standardise TCR and MHC gene symbols to their officially accepted versions, as defined by `IMGT <https://www.imgt.org/>`_.

Installation
------------

Via `PyPI <https://pypi.org/project/tidytcells/>`_ (recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``tidytcells`` can be installed using ``pip``:

.. code-block:: bash

    $ pip install tidytcells

From `source <https://github.com/yutanagano/tidytcells>`_
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The source code for the package is available `on Github <https://github.com/yutanagano/tidytcells>`_.
To install from source, clone the git repository, and run:

.. code-block:: bash

    $ pip install .

from inside the project root directory.

Useful links
------------

- `Documentation <https://tidytcells.readthedocs.io>`_
- `PyPI page <https://pypi.org/project/tidytcells>`_