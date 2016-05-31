
Command line interface
======================

The tools that can be applied to metabolic models are run through the
``psamm-model`` program. To see the full help text of the program use

.. code-block:: shell

    $ psamm-model --help

This program allows you to specify a metabolic model and a command to apply to
the given model. The available commands can be seen using the help command
given above, and are also described in more details below.

To run the program with a model, use

.. code-block:: shell

    $ psamm-model --model model.yaml command [...]

In most cases you will probably be running the command from the same directory
as where the ``model.yaml`` file is located, and in that case you can simply
run

.. code-block:: shell

    $ psamm-model command [...]

To see the help text of a command use

.. code-block:: shell

    $ psamm-model command --help

Linear programming solver
-------------------------

Many of the commands described below use a linear programming (LP) solver in
order to perform the analysis. These commands all take an option `--solver`
which can be used to select which solver to use and to specify additional
options for the LP solver. For example, in order to run the ``fba`` command
with the QSopt_ex solver, the option ``--solver name=qsoptex`` can be added:

.. code-block:: shell

    $ psamm-model fba --solver name=qsoptex

The ``--solver`` option can also be used to specify additional options for the
solver in use. For example, the Cplex solver recognizes the ``threads``
option which can be used to adjust the maximum number of threads that Cplex
will use internally (by default, Cplex will use as many threads as there are
cores on the computer):

.. code-block:: shell

    $ psamm-model fba --solver threads=4

.. _commands-fba:

Flux balance analysis (``fba``)
-------------------------------

This command will try to maximize the flux of the biomass reaction defined in
the model. It is also possible to provide a different reaction on the command
line to maximize.

To run FBA use:

.. code-block:: shell

    $ psamm-model fba

or with a specific reaction:

.. code-block:: shell

    $ psamm-model fba ATPM

By default, this performs a standard FBA and the result is output as
tab-separated values with the reaction ID, the reaction flux and the reaction
equation. If the parameter ``--loop-removal`` is given, the flux of the
internal reactions is further constrained to remove internal loops. Loop
removal is more time-consuming and under normal cicumstances the biomass
reaction flux will *not* change in response to the loop removal (only internal
reaction fluxes may change). The ``--loop-removal`` option is followed by
``none`` (no loop removal), ``tfba`` (removal using thermodynamic constraints),
or ``l1min`` (L1 minimization of the fluxes). For example, the following
command performs an FBA with thermodynamic constraints:

.. code-block:: shell

    $ psamm-model fba --loop-removal=tfba

Flux variability analysis (``fva``)
-----------------------------------

This command will find the possible flux range of each reaction when the
biomass is at the maximum value. The command will use the biomass reaction
specified in the model definition, or alternatively, a reaction can be given on
the command line.

.. code-block:: shell

    $ psamm-model fva

The output of the command will show each reaction in the model along with the
minimum and maximum possible flux values as tab-separated values. ::

    PPCK    0.0     135.266721627  [...]
    PTAr    62.3091585921    1000.0  [...]

In this example the ``PPCK`` reaction has a minimum flux of zero and maximum
flux of 135.3 units. The ``PTAr`` reaction has a minimum flux of 62.3 and a
maximum of 1000 units.

If the parameter ``--tfba`` is given, additonal thermodynamic constraints will
be imposed when evaluating model fluxes. This automatically removes internal
flux loops but is much more time-consuming.

Robustness (``robustness``)
---------------------------

Given a reaction to maximize and a reaction to vary, the robustness analysis
will run flux balance analysis and flux minimization while fixing the reaction
to vary at each iteration. The reaction will be fixed at a given number of
steps between the minimum and maximum flux value specified in the model.

.. code-block:: shell

    $ psamm-model robustness \
        --steps 200 --minimum -20 --maximum 160 EX_Oxygen

In the example above, the biomass reaction will be maximized while the
``EX_Oxygen`` (oxygen exchange) reaction is fixed at a certain flux in each
iteration. The fixed flux will vary between the minimum and maximum flux. The
number of iterations can be set using ``--steps``. In each iteration, all
reactions and the corresponding fluxes will be shown in a table, as well as
the value of the fixed flux. If the fixed flux results in an infeasible model,
no output will be shown for that iteration.

The output of the command is a list of tab-separated values indicating a
reaction ID, the flux of the varying reaction, and the flux of the reaction
with the given ID.

If the parameter ``--loop-removal`` is given, additional constraints on the
model can be imposed that remove internal flux loops. See the section on the
:ref:`commands-fba` command for more information on this option.

Random sparse network (``randomsparse``)
----------------------------------------

Delete reactions randomly until the flux of the biomass reaction falls below
the threshold. Keep deleting reactions until no more reactions can be deleted.
This can also be applied to other reactions than the biomass reaction by
specifying the reaction explicitly.

.. code-block:: shell

    $ psamm-model randomsparse 95%

When the given reaction is the biomass reaction, this results in a smaller
model which is still producing biomass within the tolerance given by the
threshold. The tolerance can be specified as a relative value (as above) or as
an absolute flux. Aggregating the results from multiple random sparse networks
allows classifying reactions as essential, semi-essential or non-essential.

If the option ``--exchange`` is given, the model will only try to delete
exchange reactions. This can be used to provide putative minimal media for
the model.

The output of the command is a tab-separated list of reaction IDs and a value
indicating whether the reaction was eliminated (``0`` when eliminated, ``1``
otherwise). If multiple minimal networks are desired, the command can be run
again and it will sample another random minimal network.

Gene Deletion (``genedelete``)
----------------------------------------

Delete single and multiple genes from a model. Once gene(s) are given the
command will delete reactions from the model requiring the gene(s) specified.
The reactions deleted will be returned as a set as well as the flux of the
model with the specified gene(s) removed.

.. code-block:: shell

    $ psamm-model genedelete

To delete genes the option ``--gene`` must be entered followed by the desired
gene ID specified in the model files. To delete multiple genes, each new gene
must first be followed by a ``--gene`` command. For example:

.. code-block:: shell

    $ psamm-model genedelete --gene ExGene1 --gene ExGene2

Single and multiple gene deletions can also be added as a text file to the
directory of the model. This can allow for simple calling of desired gene
deletions without repeating entries into the command line. In the text file
containing the desired gene deletion(s) ensure that there is one gene
ID per line. For example:

.. code-block:: shell

    $ psamm-model genedelete --gene @gene_file.txt

.. code-block:: shell

    gene_file.txt
    1 ExGene1
    2 ExGene2

Flux coupling analysis (``fluxcoupling``)
-----------------------------------------

The flux coupling analysis identifies any reaction pairs where the flux of one
reaction constrains the flux of another reaction. The reactions can be coupled
in three distinct ways depending on the ratio between the reaction fluxes.

The reactions can be fully coupled (the ratio is static and non-zero);
partially coupled (the ratio is bounded and non-zero); or directionally
coupled (the ratio is non-zero).

.. code-block:: shell

    $ psamm-model fluxcoupling

Stoichiometric consistency check (``masscheck``)
------------------------------------------------

A model or reaction database can be checked for stoichiometric inconsistencies
(mass inconsistencies). The basic idea is that we should be able to assign a
positive mass to each compound in the model and have each reaction be balanced
with respect to these mass assignments. If it can be shown that assigning the
masses is impossible, we have discovered an inconsistency.

Some variants of this idea is implemented in the :mod:`psamm.massconsistency`
module. The mass consistency check can be run using

.. code-block:: shell

    $ psamm-model masscheck

This will first try to assign a positive mass to as many compounds as possible.
This will indicate whether or not the model is consistent but in case it is
*not* consistent it is often hard to figure out how to fix the model from this
list of masses::

    [...]
    INFO: Checking stoichiometric consistency of reactions...
    C0223	1.0	Dihydrobiopterin
    C9779	1.0	2-hydroxy-Octadec-ACP(n-C18:1)
    EC0065	0.0	H+[e]
    C0065	0.0	H+
    INFO: Consistent compounds: 832/834

In this case the `H+` compounds were inconsistent because they were not
assigned a non-zero mass. A different check can be run where the residual mass
is minimized for all reactions in the model. This will often give a better idea
of which reactions need fixing::

.. code-block:: shell

    $ psamm-model masscheck --type=reaction

The following output might be generated from this command::

    [...]
    INFO: Checking stoichiometric consistency of reactions...
    IR01815	7.0     (6) |H+[c]| + |Uroporphyrinogen III[c]| [...]
    IR00307	1.0     |H+[c]| + |L-Arginine[c]| => [...]
    IR00146	0.5     |UTP[c]| + |D-Glucose 1-phosphate[c]| => [...]
    [...]
    INFO: Consistent reactions: 946/959

This is a list of reactions with non-zero residuals and their residual values.
In the example above the three reactions that are shown have been assigned a
non-zero residual (7, 1 and 0.5, respectively). This means that there is an
issue either with this reaction itself or a closely related one. In this
example the first two reactions were missing a number of `H+` compounds for
the reaction to balance.

Now the mass check can be run again marking the reactions above as checked::

    $ psamm-model masscheck --type=reaction --checked IR01815 \
        --checked IR00307 --checked IR00146
    [...]
    IR00149 0.5     |ATP[c]| + |D-Glucose[c]| => [...]

The output has now changed and the remaining residual has been shifted to
another reaction. This iterative procedure can be continued until all
stoichiometric inconsistencies have been corrected. In this example the
`IR00149` reaction also had a missing `H+` for the reaction to balance. After
fixing this error the model is consistent and the `H+` compounds can be
assigned a non-zero mass::

    $ psamm-model masscheck
    [...]
    EC0065	1.0	H+[e]
    C0065	1.0	H+
    INFO: Consistent compounds: 834/834

Formula consistency check (``formulacheck``)
--------------------------------------------

Similarly, a model or reaction database can be checked for formula
inconsistencies when the chemical formulae of the compounds in the model are
known.

.. code-block:: shell

    $ psamm-model formulacheck

For each inconsistent reaction, the reaction identifier will be printed
followed by the elements ("atoms") in, respectively, the left- and right-hand
side of the reaction, followed by the elements needed to balance the left- and
right-hand side, respectively.

Charge consistency check (``chargecheck``)
------------------------------------------

The charge check will evaluate whether the compound charge is balanced in all
reactions of the model. Any reactions that have an imbalance of charge will be
reported along with the excess charge.

.. code-block:: shell

    $ psamm-model chargecheck

Flux consistency check (``fluxcheck``)
--------------------------------------

The flux consistency check will report any reactions that are unable to take on
a non-zero flux. This is useful for finding any reactions that do not
contribute anything to the model simulation. This may indicate that the
reaction is part of a pathway that is incompletely modeled.

.. code-block:: shell

    $ psamm-model fluxcheck

If the parameter ``--tfba`` is given, additional thermodynamic constraints are
imposed when considering whether reactions can take a non-zero flux. This
automatically removes internal flux loops but is also much more time-consuming.

GapFind/GapFill (``gapfill``)
-----------------------------

The GapFind algorithm can be used to identify the compounds that are needed by
reactions in the model but cannot be produced in the model. The GapFill
algorithm will try to compute an extension of the model with reactions from the
reaction database and try to find a minimal subset that allows all blocked
compounds to be produced. This command will run GapFind to identify the blocked
compounds and then uses GapFill to try to reconstruct a model that allows these
compounds to be produced.

.. code-block:: shell

    $ psamm-model gapfill

The command will first output a list of blocked compounds and then it will list
the suggested reactions to add the model in order to unblock the blocked
compounds.

These algorithms are defined in terms of MILP problems and are therefore
(particularly GapFill) computationally expensive to run for larger models.

FastGapFill (``fastgapfill``)
-----------------------------

The FastGapFill algorithm tries to reconstruct a flux consistent model (i.e. a
model where every reaction takes a non-zero flux for at least one solutions).
This is done by extending the model with reactions from the reaction database
and trying to find a minimal subset that is flux consistent. The solution is
approximate.

The database reactions can be assigned a weight (or "cost") using the
``--penalty`` option. These weights are taken into account when determining the
minimal solution.

.. code-block:: shell

    $ psamm-model fastgapfill --penalty penalty.tsv

SBML Export (``sbmlexport``)
----------------------------

Exports the model to the SBML file format.

.. code-block:: shell

    $ psamm-model sbmlexport > model.xml

Search (``search``)
-------------------

This command can be used to search in a database for compounds or reactions. To
search for a compound use

.. code-block:: shell

    $ psamm-model search compound [...]

Use the ``--name`` option to search for a compound with a specific name or use
the ``--id`` option to search for a compound with a specific identifier.

To search for a reaction use

.. code-block:: shell

    $ psamm-model search reaction [...]

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

    $ psamm-model console
