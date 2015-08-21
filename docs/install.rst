
Install
=======

The Python module can be installed using ``pip``. This will typically require
*root* permissions.

.. code-block:: shell

    $ pip install psamm

Another option that does not require *root* permissions is to use a
`Virtualenv`_. First set up a new environment in your project directory and
activate it:

.. code-block:: shell

    $ virtualenv env
    $ . env/bin/activate

Now the Python module can be installed in the virtual environment using the
``pip`` command without requiring *root* permissions. When returning to the
project, simply reactivate the environment by running the second command.

The *psamm-import* tool is developed in a separate Git repository. After
installing PSAMM, the *psamm-import* tool can be installed using:

.. code-block:: shell

    $ pip install git+https://github.com/zhanglab/psamm-import.git

Dependencies
------------

- Linear programming solver (*Cplex*, *Gurobi*, or *QSopt_ex*)
- PyYAML (for reading the native model format)
- NumPy (optional; model matrix can be exported to NumPy matrix if available)

PyYAML is installed automatically when PSAMM is installed through ``pip``. The
linear programming solver is not strictly required but most analyses require
one to work. The LP solver *Cplex* is the preferred solver. We recently added
support for the LP solver *Gurobi*.

The rational solver *QSopt_ex* does not support MILP problems which means that
some analyses require one of the other solvers.

Cplex
-----

The Cplex Python bindings will have to be installed manually. If you are using
a virtual environment (as described above) this should be done after activating
the virtual environment:

1. Locate the directory where Cplex was installed (e.g. ``/path/to/IBM/ILOG/CPLEX_StudioXXX``).
2. Locate the appropriate subdirectory based on your platform:
   ``cplex/python/<platform>`` (e.g. ``cplex/python/x86-64_osx``).
3. Use ``pip`` to install the package from this directory: ``pip install /path/to/IBM/ILOG/CPLEX_StudioXXX/cplex/python/<platform>``

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

QSopt_ex
--------

QSopt_ex is supported through `python-qsoptex`_ which requires `GnuMP`_ and
the `QSopt_ex library`_. After installing these libraries the Python bindings
can be installed using ``pip``:

.. code-block:: shell

    $ pip install python-qsoptex

.. _Virtualenv: https://virtualenv.pypa.io/
.. _python-qsoptex: https://pypi.python.org/pypi/python-qsoptex
.. _GnuMP: https://gmplib.org/
.. _QSopt_ex library: https://github.com/jonls/qsopt-ex
