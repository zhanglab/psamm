
Installation and Materials
==========================

This tutorial will show you how to get PSAMM up and running on your computer,
how to work with the PSAMM YAML format, how to import published models into
PSAMM, and how to apply the main tools included with PSAMM to your models.

.. contents::
   :depth: 1
   :local:

Downloading the PSAMM Tutorial Data
___________________________________

The PSAMM tutorial materials are available in the psamm-tutorial GitHub repository

These files can be downloaded using the following command:

.. code-block:: shell

    $ git clone https://github.com/zhanglab/psamm-tutorial.git

This will create a directory named ``psamm-tutorial`` in your current working
folder. You can then navigate to this directory using the following command:

.. code-block:: shell

    $ cd psamm-tutorial

Now you should be in the ``psamm-tutorial`` folder and should see the following
folders:

.. code-block:: shell

    additional_files/
    E_coli_sbml/
    E_coli_excel/
    E_coli_json/

These directories include all of the files that will be needed to run the tutorial.

PSAMM Installation
__________________

PSAMM can be installed using the Python package installer pip. We recommend
that all installations be performed under a virtual Python environment. Major
programs and dependencies include: ``psamm-model``, which supports model
checking, model simulation, and model exports; Linear programming (LP) solvers
(e.g. CPLEX, Gurobi, QSopt_ex), which provide the solution of linear
programming problems; ``psamm-import``, which supports the import of models
from SBML, JSON, and Excel formats.

Setting up a Virtual Python Environment (Virtualenv)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is recommended that the PSAMM software and dependencies should be
installed under a virtual Python environment. This can be done by using
the Virtualenv_ software. Virtualenv will set up a Python environment that
permits you to install Python packages in a local directory that will not
interfere with other programs in the global Python. The virtual environment
can be set up at any local directory that you have write permission to. For
example, here we will set up the virtual environment under the main directory
of this PSAMM tutorial. First, run the following command if you are not in
the ``psamm-tutorial`` folder:

.. _Virtualenv: https://virtualenv.pypa.io/

.. code-block:: shell

    $ cd <PATH>/psamm-tutorial

In this command, ``<PATH>`` should be substituted by the directory path to
where you created the ``psamm-tutorial``. This will change your current
directory to the ``psamm-tutorial`` directory. Then, you can create a virtual
environment in the ``psamm-tutorial`` directory:

.. code-block:: shell

    $ python3 -m venv psamm-env

That will set up the virtual environment in a folder called ``psamm-env/``.
The next step is to activate the virtual environment so that the Python that is
being used will be the one that is in the virtualenv. To do this use the
following command:

.. code-block:: shell

    $ source psamm-env/bin/activate

This will change your command prompt to the following:

.. code-block:: shell

    (psamm-env) $

This indicates that the virtual environment is activated, and any installation
of Python packages will now be installed in the virtual environment. It is
important to note that when you leave the environment and return at a later
time, you will have to reactivate the environment (use the ``source`` command
above) to be able to use any packages installed in it.

.. note::

    For Windows users, the virtual environment is installed in a different
    file structure. The ``activate`` script on these systems will reside in a
    ``Scripts`` folder. To activate the environment on these systems use the
    command:

    .. code-block:: batch

        > psamm-env\Scripts\activate

.. note::

    After activating the environment, the command ``pip list`` can be used to
    quickly get an overview of the packages installed in the environment and
    the version of each package.

Setting up a Virtual Python Environment (Anaconda)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Anaconda_ is an open-source program that allows you to create virtual
environments and download Python packages. Unlike VirtualEnv, which is a
environment manager for Python, Anaconda is both a package and an environment
manager for any programming language. Anaconda manages a list of environments
for you, making it easy to work with. Instructions on how to install Anaconda
can be found `here <https://docs.anaconda.com/anaconda/install/>`_.

.. _Anaconda: https://www.anaconda.com

To create a conda environment, you do not have to be in the ``psamm-tutorial``
directory. You can create the environment from anywhere in your system with a
specific version of Python, even if it is not pre-installed:

.. code-block:: shell

    $ conda create --name psamm-env python=<version>

Unlike VirtualEnv, there will be no ``psamm-env/`` folder. A conda environment
is not dependent on your current working directory and can be activated from
anywhere using the command:

.. code-block:: shell

    $ conda activate psamm-env

When you leave the environment and return at a later time, you will
have to reactivate the environment (use the ``conda activate`` command
above) to be able to use any packages installed in it.

.. note::

    After activating the environment, the command ``conda list`` can be used to
    quickly get an overview of the packages installed in the environment and
    the version of each package.


Installation of ``psamm-model`` and ``psamm-import``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The next step will be to install ``psamm-model`` and ``psamm-import`` as well
as their requirements. To do this, you can use the Python Package Installer,
`pip`. To install both ``psamm-import`` and ``psamm-model`` you can use the
following command:

.. code-block:: shell

    (psamm-env) $ pip install git+https://github.com/zhanglab/psamm-import.git

This will install ``psamm-import`` from its Git repository and also install its
Python dependencies automatically. One of these dependencies is
``psamm-model``, so when ``psamm-import`` is installed you will also be
installing ``psamm-model``.

If you only want to install ``psamm-model`` in the environment you can run
the following command:

.. code-block:: shell

    (psamm-env) $ pip install psamm

It is important to note that if only ``psamm-model`` is installed you will be
able to apply PSAMM only on models that are represented in the YAML format
which will be described later on in the tutorial.

Installation of LP Solvers
~~~~~~~~~~~~~~~~~~~~~~~~~~

The LP (linear programming) solvers are necessary for analysis of metabolic
fluxes using the constraint-based modeling approaches.

CPLEX is the recommended solver for PSAMM and is available with an academic
license from IBM. Make sure that you use at least **CPLEX version 12.6**.
Instructions on how to install CPLEX can be found
`here <http://www-01.ibm.com/support/docview.wss?uid=swg21444285>`_.

Once CPLEX is installed, you need to install the Python bindings under
the psamm-env virtual environment using the following command:

.. code-block:: shell

    (psamm-env) $ pip install <PATH>/IBM/ILOG/CPLEX_Studio<XXX>/cplex/python/<python_version>/<platform>


The directory path in the above command should be replaced with the path to
the IBM CPLEX installation in your computer. This will install the Python
bindings for CPLEX into the virtual environment.

.. note::
    While the CPLEX software will be installed globally, the Python bindings
    should be installed specifically under the virtual environment with the
    PSAMM installation.

PSAMM also supports the use of three other linear programming solvers, Gurobi,
QSopt_ex, and GLPK. To install the Gurobi solver, Gurobi will first need to be installed
on your computer. Gurobi can be obtained with an academic license from
here: `Gurobi`_

Once Gurobi is installed the Python bindings will need to be installed in the
virtual environment by running the setup.py script in the package directory. An
example of how this could be done on a macOS is (on other platforms the path
will be different):

.. _Gurobi: http://www.gurobi.com/registration/download-reg

.. code-block:: shell

    (psamm-env) $ cd /Library/gurobi604/mac64/
    (psamm-env) $ python setup.py install

The QSopt_ex solver can also be used with PSAMM. To install this solver, Python 3.5 or
higher is required. You will first need to install Qsopt_ex on your computer and
afterwards the Python bindings (`python-qsoptex`) can be installed in the virtual
environment:

.. code-block:: shell

    (env) $ pip install cython
    (env) $ pip install python-qsoptex

Please see the `python-qsoptex documentation`_ for more information on
installing both the library and the Python bindings.

.. _python-qsoptex_higherPython: https://github.com/jonls/python-qsoptex.git
.. _QSopt_ex_higherPython: https://github.com/jonls/qsopt-ex.git

.. note::
    The QSopt_ex solver does not support Integer LP problems and as a result
    cannot be used to perform flux analysis with thermodynamic constraints. If this
    solver is used thermodynamic constraints cannot be used during simulation. By
    default ``psamm-model`` will not use these constraints.

The GLPK solver is also supported by PSAMM. The GLPK library can be installed in
the virtual environment using the following command:

.. code-block:: shell

    (psamm-env) $ pip install swiglpk

Once a solver is installed you should now be able to fully use all of the
``psamm-model`` flux analysis functions. To see a list of the installed solvers
the use the ``psamm-list-lpsolvers`` command:

.. code-block:: shell

    (psamm-env) $ psamm-list-lpsolvers

You will see the details on what solvers are installed currently and
avaliable to PSAMM. For example if the Gurobi and CPLEX solvers were both
installed you would see the following output from ``psamm-list-lpsolvers``:

.. code-block:: shell

    Prioritized solvers:
    Name: cplex
    Priority: 10
    MILP (integer) problem support: True
    QP (quadratic) problem support: True
    Rational solution: False
    Class: <class 'psamm.lpsolver.cplex.Solver'>

    Name: gurobi
    Priority: 9
    MILP (integer) problem support: True
    QP (quadratic) problem support: False
    Rational solution: False
    Class: <class 'psamm.lpsolver.gurobi.Solver'>

    Name: glpk
    Priority: 8
    MILP (integer) problem support: True
    QP (quadratic) problem support: False
    Rational solution: False
    Class: <class 'psamm.lpsolver.glpk.Solver'>

    Unavailable solvers:
    qsoptex: Error loading solver: No module named 'qsoptex'

By default the solver with the highest priority (highest priority number) is
used in constraint based simulations. If you want to use a solver with a
lower priority you will need to specify it by using the ``--solver`` option.
For example to run FBA on a model while using the Gurobi solver the following
command would be used:

.. code-block:: shell

    (psamm-env) $ psamm-model fba --solver name=gurobi

PSAMM Model Collection
______________________

Converted versions of 57 published SBML metabolic models, 9 published Excel
models and one MATLAB model are available in the `PSAMM Model Collection`_ on
GitHub. These models were converted to the YAML format and then manually edited
when needed to produce models that can generate non-zero biomass fluxes. The
changes to the models are tracked in the Git history of the repository so you
can see exactly what changes needed to be made to the models. To download
and use these models with `PSAMM` you can clone the Git repository using the
following command:

.. _`PSAMM Model Collection`: https://github.com/zhanglab/psamm-model-collection

.. code-block:: shell

    $ git clone https://github.com/zhanglab/psamm-model-collection.git

This will create the directory ``psamm-model-collection`` in your current
folder that contains one directory named ``excel`` with the converted Excel
models, one directory named ``sbml`` with the converted SBML models and one
named ``matlab`` with the converted MATLAB model. These models can then be
used for simulations with `PSAMM` using the commands detailed in this tutorial.