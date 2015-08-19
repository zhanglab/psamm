
FAQ
===

**When I run PSAMM it exits with the error "No solvers available". How can I
fix this?**

This means that PSAMM is searching for a linear programming solver but was not
able to find one. This can occur even when the Cplex solver is installed
because the Cplex Python-bindings have to be installed separately from the
main Cplex package (see :ref:`install-cplex`). Also, if using a `Virtualenv`
with Python, the Cplex Python-bindings must be installed `in the Virtualenv`.
The bindings will `not` be available in the Virtualenv if they are installed
globally.

To check whether the Cplex Python-bindings are correctly installed use the
following command:

.. code-block:: shell

    (env) $ pip list

This will show a list of Python packages that are available. The package
``cplex`` will appear in this list if the Cplex Python-bindings are correctly
installed. If the package does `not` appear in the list, follow the instuctions
at :ref:`install-cplex` to install the package.
