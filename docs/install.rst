
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

The *psamm-import* tool is developed in a separate Git repository. After
installing PSAMM, the *psamm-import* tool can be installed using:

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
   (e.g. ``cplex/python/2.7/x86-64_osx``).
3. Use ``pip`` to install the package from this directory using the following
   command.

.. code-block:: shell

    (env) $ pip install \
        /path/to/IBM/ILOG/CPLEX_Studio1262/cplex/python/<version>/<platform>

Further documentation on installing Cplex can be found in
`the Cplex documentation <http://www-01.ibm.com/support/docview.wss?uid=swg21444285>`_.

.. note::

    Python 3 support was added in a recent release of Cplex. Older versions
    only support Python 2. If you are using Python 3 make sure that you have
    the latest version of Cplex installed.

Gurobi
------

The Gurobi Python bindings will have to be installed into the virtualenv. After
activating the virtualenv:

1. Locate the directory where the Guropi python bindings were installed. For
   example, on OSX this directory is ``/Library/gurobiXXX/mac64`` where ``XXX``
   is a version code.
2. Use ``pip`` to install the package from this directory. For example:

.. code-block:: shell

    (env) $ pip install /Library/gurobi604/mac64

GLPK
----

The GLPK solver requires the GLPK library to be installed. The ``swiglpk``
Python bindings are required for PSAMM to use the GLPK library.

.. code-block:: shell

    (env) $ pip install swiglpk

QSopt_ex
--------

QSopt_ex is supported through `python-qsoptex`_ which requires `GnuMP`_ and
the `QSopt_ex library`_. After installing these libraries the Python bindings
can be installed using ``pip``:

.. code-block:: shell

    (env) $ pip install python-qsoptex

.. _Virtualenv: https://virtualenv.pypa.io/
.. _python-qsoptex: https://pypi.python.org/pypi/python-qsoptex
.. _GnuMP: https://gmplib.org/
.. _QSopt_ex library: https://github.com/jonls/qsopt-ex
