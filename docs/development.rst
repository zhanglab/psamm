
Development
===========

Test suite
----------

The python modules have test suites that allows us to automatically test
various aspects of the module implementation. It is important to make sure that
all tests run without failure *before* committing changes to any of the
modules. The test suite is run by changing to the project directory and running

.. code-block:: shell

    $ pip install -e .
    $ pip install pytest
    $ pytest

Tests that require an LP solver are skipped when no solver is installed. The
environment variable ``PSAMM_SOLVER`` can be used to select which solver is
used for the tests, e.g. to run the tests with GLPK:

.. code-block:: shell

    $ PSAMM_SOLVER=glpk pytest

Adding new tests
----------------

Adding or improving tests for python modules is highly encouraged. A test suite
for a new module should be created in ``tests/test_<modulename>.py``. These
test suites use the built-in :mod:`unittest` module.

Documentation tests
-------------------

In addition, some modules have documentation that can be tested using the
:mod:`doctest` module. These test suites should also run without failure
before any commits. They can be run by specifying the particular module (e.g
the ``affine`` module in ``expression``) using

.. code-block:: shell

    $ python -m psamm.expression.affine -v
