
Command line interface
======================

The tools that can be applied to metabolic models are run through the
``model`` program. To see the full help text of the program use

.. code-block:: shell

    $ model --help

This program allows you to specify a metabolic model and a command to apply to
the given model. The available commands can be seen using the help command
given above, and are also described in more details below.

To run the program with a model, use

.. code-block:: shell

    $ model --model model.yaml command [...]

In most cases you will probably be running the command from the same directory
as where the ``model.yaml`` file is located, and in that case you can simply
run

.. code-block:: shell

    $ model command [...]

To see the help text of a command use

.. code-block:: shell

    $ model command --help

Flux balance analysis (``fba``)
-------------------------------

This command will first try to maximize the flux of the biomass reaction
defined in the model. It is also possible to provide a different reaction on
the command line to maximize.

By default, this is followed by running the flux balance analysis with
thermodynamic constraints (tFBA) in order to remove internal flux cycles. The
results of these two analyses is presented in a two-column table along with the
reaction IDs.

If the parameter ``--no-tfba`` is given, the second column instead represents a
flux minimization in which the FBA maximum is fixed while the sum of the fluxes
is minimized. This will often eliminate loops as well.

To run FBA use:

.. code-block:: shell

    $ model fba

or with a specific reaction:

.. code-block:: shell

    $ model fba ATPM

Robustness (``robustness``)
---------------------------

Given a reaction to maximize and a reaction to vary, the robustness analysis
will run flux balance analysis and flux minimization while fixing the reaction
to vary at each iteration. The reaction will be fixed at a given number of
steps between the minimum and maximum flux value specified in the model.

.. code-block:: shell

    $ model robustness \
	   --steps 200 --minimum -20 --maximum 160 EX_Oxygen

In the example above, the biomass reaction will be maximized while the
``EX_Oxygen`` (oxygen exchange) reaction is fixed at a certain flux in each
iteration. The fixed flux will vary between the minimum and maximum flux. The
number of iterations can be set using ``--steps``. In each iteration, all
reactions and the corresponding fluxes will be shown in a table, as well as
the value of the fixed flux. If the fixed flux results in an infeasible model,
no output will be shown for that iteration.

Random sparse network (``randomsparse``)
----------------------------------------

Delete reactions randomly until the flux of the biomass reaction falls below
the threshold. Keep deleting reactions until no more reactions can be deleted.
This can also be applied to other reactions than the biomass reaction by
specifying the reaction explicitly.

.. code-block:: shell

    $ model randomsparse 0.95

When the given reaction is the biomass reaction, this results in a smaller
model which is still producing biomass within the tolerance given by the
threshold. Aggregating the results from multiple random sparse networks allows
classifying reactions as essential, semi-essential or non-essential.

Mass consistency check (``masscheck``)
--------------------------------------

A model or reaction database can be checked for mass inconsistencies. The basic
idea is that we should be able to assign a positive mass to each compound in the
model and have each reaction be balanced with respect to these mass assignments.
If it can be shown that assigning the masses is impossible, we have discovered
an inconsistency.

Some variants of this idea is implemented in the :mod:`metnet.massconsistency`
module. The mass consistency check can be run using

.. code-block:: shell

    $ model masscheck

This will first try to assign a positive mass to as many compounds as possible.
This will indicate whether or not the model is consistent but in case it is
*not* consistent it is often hard to figure out how to fix the model from this
list of masses.

Next, a different check is run where the residual mass is minimized for all
reactions in the model. This will often give a better idea of which reactions
need fixing.

Formula consistency check (``formulacheck``)
--------------------------------------------

Similarly, a model or reaction database can be checked for formula
inconsistencies when the chemical formulae of the compounds in the model are
known.

.. code-block:: shell

    $ model formulacheck

For each inconsistent reaction, the reaction identifier will be printed
followed by the elements ("atoms") in, respectively, the left- and right-hand
side of the reaction, followed by the elements needed to balance the left- and
right-hand side, respectively.

GapFind/GapFill (``gapfill``)
-----------------------------

The GapFind algorithms can be used to identify the compounds that are needed by
reactions in the model but cannot be produced in the model. The GapFill
algorithm will extend the model with reactions from the parent database and try
to find a minimal subset that allows all blocked compounds to be produced. This
command will run GapFind to identify the blocked compounds and then uses
GapFill to try to reconstruct a model that allows these compounds to be
produced.

These algorithms are defined in terms of MILP problems and are therefore
(particularly GapFill) computationally expensive to run for larger models.

.. code-block:: shell

    $ model gapfill

FastGapFill (``fastgapfill``)
-----------------------------

The FastGapFill algorithm tries to reconstruct a flux consistent model (i.e. a
model where every reaction takes a non-zero flux for at least one solutions).
This is done by extending the model with reactions from the parent database and
trying to find a minimal subset that is flux consistent. The solution is
approximate.

The database reactions can be assigned a weight (or "cost") using the
``--penalty`` option. These weights are taken into account when determining the
minimal solution.

.. code-block:: shell

    $ model fastgapfill --penalty penalty.tsv

Search (``search``)
-------------------

This command can be used to search in a database for compounds or reactions. To
search for a compound use

.. code-block:: shell

    $ model search compound [...]

Use the ``--name`` option to search for a compound with a specific name or use
the ``--id`` option to search for a compound with a specific identifier.

To search for a reaction use

.. code-block:: shell

    $ model search reaction [...]

Use the ``--id`` option to search for a reaction with a specific identifier.
The ``--compound`` option can be used to search for reactions that include a
specific compound. If more that one compound identifier is given
(comma-separated) this will find reactions that include all of the given
compounds.

Console (``console``)
---------------------

This command will start a Python session where the model has been loaded into
the corresponding Python object representation.

.. code-block:: shell

    $ model console
