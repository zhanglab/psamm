
Install
=======

PSAMM can be installed using the Python package installer ``pip``. We recommend
that you use a `Virtualenv`_ when installing PSAMM. First, create a Virtualenv
in your project directory and activate the environment. On Linux/OSX the
following terminal commands can be used:

.. code-block:: shell

    $ virtualenv env
    $ source env/bin/activate

Then, install PSAMM using the ``pip`` command:

.. code-block:: shell

    (env) $ pip install psamm

When returning to the project from a new terminal window, simply reactivate
the environment by running

.. code-block:: shell

    $ source env/bin/activate

The *psamm-import* tool is included in the main PSAMM repository. Some
additional model specific importers for Excel format models associated
with publications are maintained in a separate repository. After
installing PSAMM, support for these import functions can be added through
installing this additional program:

.. code-block:: shell

    (env) $ pip install git+https://github.com/zhanglab/psamm-import.git

Dependencies
------------

- Linear programming solver (*Cplex*, *Gurobi*, *GLPK* or *QSopt_ex*)
- PyYAML (for reading the native model format)
- NumPy (optional; model matrix can be exported to NumPy matrix if available)

PyYAML is installed automatically when PSAMM is installed through ``pip``. The
linear programming solver is not strictly required but most analyses require
one to work. The LP solver *Cplex* is the preferred solver. We recently added
support for the LP solver *Gurobi* and *GLPK*.

The rational solver *QSopt_ex* does not support MILP problems which means that
some analyses require one of the other solvers. The MILP support in *GLPK* is
still experimental so it is disabled by default.

.. _install-cplex:

Cplex
-----

The Cplex Python bindings will have to be installed manually. Make sure that
you are using at least **Cplex version 12.6**. If you are using
a virtual environment (as described above) this should be done after activating
the virtual environment:

1. Locate the directory where Cplex was installed (e.g. ``/path/to/IBM/ILOG/CPLEX_StudioXXX``).
2. Locate the appropriate subdirectory based on your platform and Python
   version: ``cplex/python/<version>/<platform>``
   (e.g. ``cplex/python/3.7/x86-64_osx``).
3. Use ``pip`` to install the package from this directory using the following
   command.

.. code-block:: shell

    (env) $ pip install \
        /path/to/IBM/ILOG/CPLEX_Studio1262/cplex/python/<version>/<platform>

Further documentation on installing Cplex can be found in
`the Cplex documentation <http://www-01.ibm.com/support/docview.wss?uid=swg21444285>`_.


Gurobi
------

The Gurobi Python bindings will have to be installed into the virtualenv. After
activating the virtualenv:

1. Change the directory to where the Guropi python bindings were installed. For
   example, on OSX this directory is ``/Library/gurobiXXX/mac64`` where ``XXX``
   is a version code.
2. Use ``python`` to install the package from this directory. For example:

.. code-block:: shell

    (env) $ cd /Library/gurobi604/mac64
    (env) $ python setup.py install

GLPK
----

The GLPK solver requires the GLPK library to be installed. The ``swiglpk``
Python bindings are required for PSAMM to use the GLPK library.

.. code-block:: shell

    (env) $ pip install swiglpk

QSopt_ex
--------

QSopt_ex is supported through `python-qsoptex`_ which requires `python-qsoptex_higherPython`_
and `QSopt_ex_higherPython`_ . This can be installed using ``pip``:

.. code-block:: shell

    (env) $ pip install cython
    (env) $ pip install python-qsoptex

.. _python-qsoptex_higherPython: https://github.com/jonls/python-qsoptex.git
.. _QSopt_ex_higherPython: https://github.com/jonls/qsopt-ex.git

LP Solver Compatibility
-----------------------

Not all of the LP solvers supported are supported across all python versions.
A table showing which solvers are compatible with which versions of python
is shown below:

.. list-table:: Python Version Compatibility
   :header-rows: 1

   * - Solver
     - Python 3.5
     - Python 3.6
     - Python 3.7
     - Python 3.8
     - Python 3.9
   * - Cplex
     - Yes
     - Yes
     - Yes
     - Yes
     - Yes
   * - Qsopt_ex
     - Yes
     - Yes
     - Yes
     - Yes
     - Yes
   * - Gurobi
     - No
     - No
     - Yes
     - Yes
     - Yes
   * - GLPK
     - No
     - Yes
     - Yes
     - Yes
     - Yes

LP Solver Global Options
------------------------

Additionally, not all LP solvers are compatible with all of the
global parameters for solvers. Reference the table below for the
global parameters you might require before choosing a solver.


.. list-table:: Global Solver Options
   :header-rows: 1

   * - Solver
     - feasibility tolerance
     - optimality tolerance
     - integrality tolerance
     - threads
   * - Cplex
     - Yes
     - Yes
     - Yes
     - Yes
   * - Qsopt_ex
     - Yes
     - Yes
     - No
     - No
   * - Gurobi
     - Yes
     - Yes
     - Yes
     - Yes
   * - GLPK
     - Yes
     - Yes
     - Yes
     - No
