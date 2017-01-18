
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
line to maximize. [Orth10]_ [Fell86]_

To run FBA use:

.. code-block:: shell

    $ psamm-model fba

or with a specific reaction:

.. code-block:: shell

    $ psamm-model fba --objective=ATPM

By default, this performs a standard FBA and the result is output as
tab-separated values with the reaction ID, the reaction flux and the reaction
equation. If the parameter ``--loop-removal`` is given, the flux of the
internal reactions is further constrained to remove internal loops
[Schilling00]_. Loop removal is more time-consuming and under normal
circumstances the biomass reaction flux will *not* change in response to the
loop removal (only internal reaction fluxes may change). The ``--loop-removal``
option is followed by ``none`` (no loop removal), ``tfba`` (removal using
thermodynamic constraints), or ``l1min`` (L1 minimization of the fluxes). For
example, the following command performs an FBA with thermodynamic constraints:

.. code-block:: shell

    $ psamm-model fba --loop-removal=tfba

Flux variability analysis (``fva``)
-----------------------------------

This command will find the possible flux range of each reaction when the
biomass is at the maximum value [Mahadevan03]_. The command will use the
biomass reaction specified in the model definition, or alternatively, a
reaction can be given on the command line following the ``--objective`` option.

.. code-block:: shell

    $ psamm-model fva

The output of the command will show each reaction in the model along with the
minimum and maximum possible flux values as tab-separated values. ::

    PPCK    0.0     135.266721627  [...]
    PTAr    62.3091585921    1000.0  [...]

In this example the ``PPCK`` reaction has a minimum flux of zero and maximum
flux of 135.3 units. The ``PTAr`` reaction has a minimum flux of 62.3 and a
maximum of 1000 units.

If the parameter ``--loop-removal=tfba`` is given, additonal thermodynamic
constraints will be imposed when evaluating model fluxes. This automatically
removes internal flux loops [Schilling00]_ but is much more time-consuming.

Robustness (``robustness``)
---------------------------

Given a reaction to maximize and a reaction to vary, the robustness analysis
will run flux balance analysis and flux minimization while fixing the reaction
to vary at each iteration. The reaction will be fixed at a given number of
steps between the minimum and maximum flux value specified in the model
[Edwards00]_.

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
must first be followed by a ``--gene`` option. For example:

.. code-block:: shell

    $ psamm-model genedelete --gene ExGene1 --gene ExGene2

The list of genes to delete can also be specified in a text file. This allows
to you perform many gene deletions by simply specifying the file name when
running the ``genedelete`` command. The text file must contain one gene ID per
line. For example:

.. code-block:: shell

    $ psamm-model genedelete --gene @gene_file.txt

The file gene_file.txt would contain the following lines::

    ExGene1
    ExGene2

To delete genes using different algorithms use ``--method`` to specify
which algorithm for the solver to use. The default implementation for this
command is FBA. The two examples below will yield the same result.

.. code-block:: shell

    $ psamm-model genedelete --gene ExGene1

.. code-block:: shell

    $ psamm-model genedelete --gene ExGene1 --method fba

To delete genes using the Minimization of Metabolic Adjustment (MOMA)
algorithm use the command line argument ``--method moma``.

.. code-block:: shell

    $ psamm-model genedelete --gene ExGene1 --method moma

There are four implementations of MOMA inside of PSAMM:

lin MOMA (``--method lin_moma``)
  Finds the maximum biomass after a gene deletions, such that the change in
  the flux system is minimized when compared to the wild type. Minimization
  is done by minimizing the \|wild type fluxes - knockout fluxes|. This helps
  avoid the assumption that an organism will perform optimally directly after
  removing a gene.

MOMA (``--method moma``) : Recommended
  Finds the maximum biomass after a gene deletions, such that the change in
  the flux system is minimized when compared to the wild type. Minimization
  is done by minimizing the (wild type fluxes - knockout fluxes) :sup:`2` .
  This helps avoid the assumption that an organism will perform optimally
  directly after removing a gene.

lin MOMA 2 (``--method lin_moma2``) : Experimental
  Similar to ``lin_moma``, but this implementation solves for the wild type
  fluxes at the same time as the knockout fluxes to ensure not to rely on the
  arbitrary flux vector found with FBA. This will find a more optimal solution
  to the problem than the original ``lin_moma``.

MOMA 2 (``--method moma2``) : Experimental
  Similar to ``moma``, but this implementation solves for the wild type
  fluxes at the same time as the knockout fluxes to ensure not to rely on the
  arbitrary flux vector found with FBA. This will find a more optimal solution
  to the problem than the original ``moma``.

Flux coupling analysis (``fluxcoupling``)
-----------------------------------------

The flux coupling analysis identifies any reaction pairs where the flux of one
reaction constrains the flux of another reaction. The reactions can be coupled
in three distinct ways depending on the ratio between the reaction fluxes
[Burgard04]_.

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
masses is impossible, we have discovered an inconsistency [Gevorgyan08]_.

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

If the parameter ``--loop-removal=tfba`` is given, additional thermodynamic
constraints are imposed when considering whether reactions can take a non-zero
flux. This automatically removes internal flux loops but is also much more
time-consuming.

Reaction duplicates check (``duplicatescheck``)
-----------------------------------------------

This command simply checks whether multiple reactions exist in the model that
have the same or similar reaction equations. By default, this check will ignore
reaction directionality and stoichiometric values when considering whether
reactions are identical. The options ``--compare-direction`` and
``--compare-stoichiometry`` can be used to make the command consider these
properties as well.

.. code-block:: shell

    $ psamm-model duplicatescheck

Gap check (``gapcheck``)
------------------------

This gap check command will try to identify the compounds in the model
that cannot be produced. This is useful for identifying incomplete pathways in
the model. The command will report a list of all compounds in the model that
are blocked for production.

.. code-block:: shell

    $ psamm-model gapcheck

When checking whether a compound can be produced, it is sufficient for
production that all precursors can be produced and it is *not* necessary for
every compound to also be consumed by another reaction (in other words, for
the purpose of this analysis there are implicit sinks for every compound in
the model). This means that even if this command reports that no compounds are
blocked, it may still not be possible for the model to be viable under the
steady-state assumption of FBA. The option ``--no-implicit-sinks`` can be used
to perform the gap check without implicit sinks.

The gap check is performed with the medium that is defined in the model. It
may be useful to run the gap check with every compound in the medium available.
This can easily be done by specifying the ``--unrestricted-exchange`` option
which removes all limits on the exchange reactions during the check.

There are some additional gap checking methods that can be enabled with the
``--method`` option. The method ``sinkcheck`` can be used to find compounds
that cannot be synthesized from scratch. The standard gap check will report
compounds as produced if they can participate in a reaction, even if the
compound itself cannot be synthesized from precursors in the medium. To find
such compounds use the ``sinkcheck``. This check will generally indicate more
compounds as blocked. Lastly, the method ``gapfind`` can be used. This method
should produce the same result as the default method but is implemented in an
alternative way that is specified in [Kumar07]_. This method is *not* used by
default because it tends to result in difficulties for the solver when used
with larger models.

GapFill (``gapfill``)
---------------------

The GapFill algorithm will try to compute an extension of the model with
reactions from the reaction database and try to find a minimal subset that
allows all blocked compounds to be produced. In addition to suggesting possible
database reactions to add to the model, the command will also suggest possible
transport and exchange reactions. The GapFill algorithm implemented in this
command is a variant of the gap-filling procedure described in [Kumar07]_.

.. code-block:: shell

    $ psamm-model gapfill

The command will first list the reactions in the model followed by the
suggested reactions to add to the model in order to unblock the blocked
compounds. If ``--allow-bounds-expansion`` is specified, the procedure may also
suggest that existing model reactions have their flux bounds widened, e.g.
making an existing irreversible reaction reversible. To unblock only specific
compounds, use the ``--compound`` option:

.. code-block:: shell

    $ psamm-model gapfill --compound leu-L[c] --compound ile-L[c]

In this example, the procedure will try to add reactions so that leucine
(``leu-L``) and isoleucine (``ile-L``) in the ``c`` compartment can be
produced. Multiple compounds can be unblocked at the same time and the list of
compounds to unblock can optionally be specified as a file by prefixing the
file name with ``@``.

.. code-block:: shell

    $ psamm-model gapfill --compound @list_of_compounds_to_unblock.tsv

The GapFind algorithm is defined in terms of a MILP problem and can therefore
be computationally expensive to run for larger models.

The original GapFill algorithm uses a solution procedure which implicitly
assumes that the model contains implicit sinks for all compounds. This means
that even with the reactions proposed by GapFill the model may need to produce
compounds that cannot be used anywhere. The implicit sinks can be disabled
with the ``--no-implicit-sinks`` option.

FastGapFill (``fastgapfill``)
-----------------------------

The FastGapFill algorithm tries to reconstruct a flux consistent model (i.e. a
model where every reaction takes a non-zero flux for at least one solutions).
This is done by extending the model with reactions from the reaction database
and trying to find a minimal subset that is flux consistent. The solution is
approximate [Thiele14]_.

The database reactions can be assigned a weight (or "cost") using the
``--penalty`` option. These weights are taken into account when determining the
minimal solution.

.. code-block:: shell

    $ psamm-model fastgapfill --penalty penalty.tsv

SBML Export (``sbmlexport``)
----------------------------

Exports the model to the SBML file format. This command exports the model as
an `SBML level 3`_ file with flux bounds, objective and gene information
encoded with `Flux Balance Constraints version 2`_.

.. code-block:: shell

    $ psamm-model sbmlexport model.xml

If the file name is omitted, the file contents will be output directly to the
screen. Using the ``--pretty`` option makes the output formatted for
readability.

.. _`SBML level 3`: http://sbml.org/Documents/Specifications
.. _`Flux Balance Constraints version 2`: http://sbml.org/Documents/Specifications/SBML_Level_3/Packages/fbc

Excel Export (``excelexport``)
------------------------------

Exports the model to the Excel file format.

.. code-block:: shell

    $ psamm-model excelexport model.xls

Table Export (``tableexport``)
------------------------------

Exports the model to the tsv file format.

.. code-block:: shell

    $ psamm-model tableexport reactions > model.tsv

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
