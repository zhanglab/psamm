
Install
=======

The Python module can be installed using ``pip``. It is also possible to
install a specific tag or branch by appending the name to the Git URL (e.g.
append ``@v0.1`` to get the tag ``v0.1``). This will typically require ``root``
permission.

.. code-block:: shell

    $ pip install git+ssh://git@github.com/zhanglab/model_script.git

Another option is to use ``virtualenv``. First set up a new environment in the
project directory, and activate it.

.. code-block:: shell

    $ virtualenv env
    $ . env/bin/activate

Now the Python module can be installed in the virtual environment using the
``pip`` command without requiring ``root`` permissions. When returning to the
project, simply reactivate the environment by running the second command.

Dependencies
------------

- Linear programming solver (*Cplex*, *QSopt_ex*)
- PyYAML (for reading the native model format)
- NumPy (optional; model matrix can be exported to NumPy matrix if available)

PyYAML is installed automatically when PSAMM is installed through ``pip``. The
linear programming solver is not strictly required but most analyses require
one to work. The LP solver *Cplex* is the preferred solver. The rational solver
*QSopt_ex* does not support MILP problems which means that some analyses
require *Cplex*.

Cplex
-----

The Cplex Python bindings will have to be installed manually. If you are using
a virtual environment (as described above) this should be done after activating
the virtual environment:

1. Go to the Cplex install directory: ``cd /path/to/Cplex``
2. Go to the appropriate subdirectory based on your platform:
   ``cd cplex/python/<platform>``
3. Run ``pip install .``

QSopt_ex
--------

QSopt_ex is supported through `python-qsoptex`_ which requires `GnuMP`_ and
the `QSopt_ex library`_. After installing these libraries the Python bindings
can be installed using ``pip``:

.. code-block:: shell

    $ pip install python-qsoptex

.. _python-qsoptex: https://pypi.python.org/pypi/python-qsoptex
.. _GnuMP: https://gmplib.org/
.. _QSopt_ex library: https://github.com/jonls/qsopt-ex
