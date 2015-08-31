
Development
===========

Test suite
----------

The python modules have test suites that allows us to automatically test
various aspects of the module implementation. It is important to make sure that
all tests run without failure *before* committing changes to any of the
modules. The test suite is run by changing to the project directory and running

.. code-block:: shell

    $ ./setup.py test

To run the tests on all the supported Python platforms with additional tests
for coding style (PEP8) and building documentation, use tox_:

.. code-block:: shell

    $ tox

When testing with tox, the local path for the Cplex python module must be
provided in the environment variables ``CPLEX_PYTHON2_PACKAGE`` and
``CPLEX_PYTHON3_PACKAGE`` for Python 2 and Python 3, respectively. For example:

.. code-block:: shell

    $ export CPLEX_PYTHON2_PACKAGE=/path/to/IBM/ILOG/CPLEX_StudioXXX/cplex/python/2.7/x86-64_osx
    $ export CPLEX_PYTHON3_PACKAGE=/path/to/IBM/ILOG/CPLEX_StudioXXX/cplex/python/3.4/x86-64_osx
    $ tox -e py27-cplex,py34-cplex

.. note::

    Python 3 support was added in a recent release of Cplex. Older versions
    only support Python 2.

Similarly, the local path to the Gurobi package must be specified in the
environment variable ``GUROBI_PYTHON_PACKAGE``:

.. code-block:: shell

    $ export GUROBI_PYTHON_PACKAGE=/Library/gurobi604/mac64
    $ tox -e py27-gurobi

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

.. _tox: https://testrun.org/tox/
