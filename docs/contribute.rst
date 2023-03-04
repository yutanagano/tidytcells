Contribute
==========

Please feel free to contribute to the project by submitting `bug reports <https://github.com/yutanagano/tidytcells/issues/new>`_ and/or `pull requests <https://github.com/yutanagano/tidytcells/pulls>`_ at ``tidytcells``' `GitHub repository <https://github.com/yutanagano/tidytcells>`_!

Setting up a development environment
------------------------------------

To set up a development environment for ``tidytcells``, please clone the git repository linked above, and then set up a virtual environment by executing the following commands:

.. code-block:: bash

    $ python -m venv .venv
    $ source .venv/bin/activate
    $ pip install -e .[dev]

This will set up a python virtual environment in which an editable copy of the ``tidytcells`` package (changes in the source code will be reflected in real time), plus the development dependencies will be installed.

Running tests
-------------

The following commands should be run from the testing environment set up as per above.

To run the complete suite of tests including the test building of the `sphinx <https://www.sphinx-doc.org/en/master/#>` docs and `pytest <https://docs.pytest.org>` tests on multiple python versions, you can use `tox <https://tox.wiki/en/latest/>`:

.. code-block:: bash

    $ tox

To execute just the pytest unit tests, run:

.. code-block:: bash

    $ python -m pytest

If you want to also see test coverage, run:

.. code-block:: bash

    $ python -m pytest --cov=src