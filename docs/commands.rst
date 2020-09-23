
Command line interface
======================

The tools that can be applied to metabolic models are run through the
``psamm-model`` program. To see the full help text of the program use

.. code-block:: shell

    $ psamm-model --help

This program allows you to specify a metabolic model and a command to apply to
the given model. The available commands can be seen using the help command
given above, and are also described in more detail below.

To run the program with a model, use the following command:

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
order to perform the analysis. These commands all take an option ``--solver``
which can be used to select which solver to use and to specify additional
options for the LP solver. For example, in order to run the ``fba`` command
with the QSopt_ex solver, the option ``--solver name=qsoptex`` can be added:

.. code-block:: shell

    $ psamm-model fba --solver name=qsoptex

The ``--solver`` option can also be used to specify additional options for the
solver in use. For example, the CPLEX solver recognizes the ``threads``
option which can be used to adjust the maximum number of threads that CPLEX
will use internally (by default, CPLEX will use as many threads as there are
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
equation. If the parameter ``--loop-removal`` is given, the fluxes of the
internal reactions are further constrained to remove internal loops
[Schilling00]_. Loop removal is more time-consuming and under normal
circumstances the biomass reaction flux will *not* change in response to the
loop removal (only internal reaction fluxes may change). The ``--loop-removal``
option is followed by ``none`` (no loop removal), ``tfba`` (removal using
thermodynamic constraints), or ``l1min`` (L1 minimization of the fluxes). For
example, the following command performs an FBA with thermodynamic constraints:

.. code-block:: shell

    $ psamm-model fba --loop-removal=tfba

By default, the output of the FBA command will only display reactions which
have non-zero fluxes. This can be overridden with the ``--all-reactions``
option to display all reactions even if they have flux values of zero.

.. code-block:: shell

    $ psamm-model fba --all-reactions

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

If the parameter ``--loop-removal=tfba`` is given, additional thermodynamic
constraints will be imposed when evaluating model fluxes. This automatically
removes internal flux loops [Schilling00]_ but is much more time-consuming.

By default, FVA is performed with the objective reaction (either the biomass
reaction or reaction given through the ``--objective`` option) fixed at its
maximum value. It is also possible allow this reaction flux to vary by a
specified amount through the ``--thershold`` option. When this option is
used the variability results will show the possible flux ranges when the
objective reaction is greater than or equal to the threshold value.

The threshold can be specified by either giving a percentage of the maximum
objective flux or by giving a defined flux value.

.. code-block:: shell

    $ psamm-model fva --threshold 90%
    or
    $ psamm-model fva --threshold 1.2

The FVA command can also be run as parallel processes to speed up simulations
done on larger models. This can be done using the ``--parallel`` option. Either
a specific number of parallel jobs can be given or 0 can be given to
automatically detect and use the maximum number of parallel processes.

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

It is also possible to print out the flux of all reactions for each step in
the robustness simulation instead of just printing the varying reaction flux.
This can be done through using the ``--all-reaction-fluxes`` option.

The Robustness command can also be run as parallel processes to speed up
simulations done on larger models. This can be done using the
``--parallel`` option. Either a specific number of parallel jobs can
be given or 0 can be given to automatically detect and use the maximum
number of parallel processes.

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

By default, the randomsparse command will perform the deletions on reactions
in the model. The ``--type`` option can be used to change this deletion to
act on the genes in the model or to act on only the set of exchange reactions.
The gene deletion option will remove a gene from a network and then assess
which reactions would be affected by that gene loss based on the provided
gene associations. The exchange reaction deletion will only delete reactions
from the set of exchange reactions in the model and can be used to generate
a putative minimal media for the model.

.. code-block:: shell

    $ psamm-model randomsparse --type genes 95%
    or
    $ psamm-model randomsparse --type exchange 95%

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
which algorithm for the solver to use. The default method for the command is
FBA. To delete genes using the Minimization of Metabolic Adjustment (MOMA)
algorithm use the command line argument ``--method moma``. MOMA is based on
the assumption that the knockout organism has not had time to adjust its gene
regulation to maximize biomass production so fluxes will be close to wildtype
fluxes.

.. code-block:: shell

    $ psamm-model genedelete --gene ExGene1 --method moma

There are four variations of MOMA available in PSAMM, defined in the following
way (where :math:`\bar{v}` is the wild type fluxes and :math:`\bar{u}` is the
knockout fluxes):

MOMA (``--method moma``)
    Finds the reaction fluxes in the knockout, such that the difference in flux
    from the wildtype is minimized. Minimization is performed with the
    Euclidean distance: :math:`\sum_j (v_j - u_j)^2`. The wildtype fluxes are
    obtained from the wildtype model (i.e. before genes are deleted) by FBA
    with L1 minimization. L1 minimization is performed on the FBA result to
    remove loops and make the result disregard internal loop fluxes. [Segre02]_

Linear MOMA (``--method lin_moma``)
    Finds the reaction fluxes in the knockout, such that the difference in flux
    from the wildtype is minimized. Minimization is performed with the
    Manhattan distance: :math:`\sum_j \|v_j - u_j\|`. The wildtype fluxes are
    obtained from the wildtype model (i.e. before genes are deleted) by FBA
    with L1 minimization. L1 minimization is performed on the FBA result to
    remove loops and make the result disregard internal loop fluxes. [Mo09]_

Combined-model MOMA (``--method moma2``) (Experimental)
    Similar to ``moma``, but this implementation solves for the wild type
    fluxes at the same time as the knockout fluxes to ensure not to rely on the
    arbitrary flux vector found with FBA.

Combined-model linear MOMA (``--method lin_moma2``) (Experimental)
    Similar to ``lin_moma``, but this implementation solves for the wild type
    fluxes at the same time as the knockout fluxes to ensure not to rely on the
    arbitrary flux vector found with FBA.

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
fixing this error, the model is consistent and the `H+` compounds can be
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

Some reactions could

Reaction duplicates check (``dupcheck``)
-----------------------------------------------

This command simply checks whether multiple reactions exist in the model that
have the same or similar reaction equations. By default, this check will ignore
reaction directionality and stoichiometric values when considering whether
reactions are identical. The options ``--compare-direction`` and
``--compare-stoichiometry`` can be used to make the command consider these
properties as well.

.. code-block:: shell

    $ psamm-model dupcheck

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

The penalty file provided should be a tab separated table that contains
reaction IDs and assigned penalty values in two columns:

.. code-block:: shell

    RXN1    10
    RXN2    15
    RXN3    20
    ....

In addition to setting penalty values directly through the ``--penalty`` argument,
default penalties for all other database reactions can be set through the
``--db-penalty`` argument. Reactions that do not have penalty values
explicitly set through the ``--penalty`` argument will be assigned
this penalty value. Similarly, penalty values can be assigned for new
exchange reactions and artificial transport reactions through the
``--ex-penalty`` and ``--tp-penalty`` arguments. An example using
all three of these arguments can be seen here:

.. code-block:: shell

    $ psamm-model fastgapfill --db-penalty 100 --tp-penalty 1000 --ex-penalty 150


Predict primary pairs (``primarypairs``)
----------------------------------------------------------------------

This command is used to predict element-transferring reactant/product pairs
in the reactions of the model. This can be used to determine the flow of
elements through reactions. Two methods for predicting the pairs are available:
`FindPrimaryPairs` (``fpp``) [Steffensen17]_ and
MapMaker (``mapmaker``) [Tervo16]_. The ``--method`` option can used to select
which prediction method to use:

.. code-block:: shell

    $ psamm-model primarypairs --method=fpp

The result is reported as a table of four columns. The first column is the
reactions ID, the second and third columns contain the compound ID of the
reactant and product. The fourth column contains the predicted transfer of
elements.

The ``primarypairs`` command will run slowly on models that contain artificial
reactions such as biomass reactions or condensed biosynthesis reactions.
Because the reactant/product pair prediction in these reactions is
not as biologically meaningful these reactions can be excluded through
the ``--exclude`` option. This option can be used to either give reaction
IDs to exclude or to give an input file containing a list of reactions
IDs to exclude:

.. code-block:: shell

    $ psamm-model primarypairs --exclude BiomassReaction
    or
    $ psamm-model primarypairs --exclude @./exclude.tsv


PSAMM-Vis (``vis``)
-----------------------------

Models can be visualized through the use of `PSAMM-vis` as implemented in the
``vis`` command in `PSAMM`. This command can use
the `FindPrimaryPairs` algorithm to help generate images of full models
or subsets of models. The output of this command will consist of a graph
file in the `dot` language, ``reactions.dot``, and two files
called ``reactions.nodes.tsv`` and ``reactions.edges.tsv`` that contain
the network data in TSV format.

To run the ``vis`` command the following command can be used:

.. code-block:: shell

    $ psamm-model vis

Basic Graph Generation
~~~~~~~~~~~~~~~~~~~~~~

By default, the ``vis`` command uses the `FindPrimaryPairs` algorithm to
simplify the graph that is produced. This algorithm runs much faster
if certain types of artificial reactions are not considered when doing
the reactant/product pair prediction. These reactions often represent
Biomass production or condensed biosynthesis processes. To exclude
these reactions the ``vis`` command can be run with the ``--exclude``
option followed by an input file that contains a list of reaction IDs:

.. code-block:: shell

    $ psamm-model vis --exlcude @{path to file}

Running this command, `PSAMM-vis` only produces three files described above.
Graph image generating softwares can convert these files to
actual images. If the program `Graphviz` is installed on the computer,
then this program can be used within `PSAMM` to generate the network image
directly. This can be done by adding the ``--image`` argument followed
by any `Graphviz` supported image format:

.. code-block:: shell

    $ psamm-model vis --image {format (pdf, eps, svg, etc.)}

In addition, biomass reaction defined in model.yaml will be excluded automatically.

While the ``vis`` function in `PSAMM` uses `FindPrimaryPairs` for graph
simplification by default, the command is also able to run using no
graph simplification (``no-fpp``).

This can be done through using the ``--method`` argument:

.. code-block:: shell

    $ psamm-model vis --method no-fpp

The resulting graphs can be further simplified to only show element
transfers that contain a specified element through the ``--element`` argument.
When using this option any reactant/product pairs that do not transfer
the specified element will not be shown on the graph. To use this option the
following command can be used:

.. code-block:: shell

    $ psamm-model vis --element {Atomic Symbol}

Additionally, the final graphs created through ``vis`` command can only show
a specified subset of reactions from the larger model. This can be done using
``--subset`` argument to provide a file containing a single column list of
either reaction IDs or metabolite IDs (but not the mix of reaction and
compound IDs).

.. code-block:: shell

    $ psamm-model vis --subset {path to file}

And example of this file would be:

.. code-block:: shell

    rxn_1
    rxn_2
    rxn_4

Or:

.. code-block:: shell

    cpd_1[c]
    cpd_1[e]
    cpd_2[c]


Further modification can be done to the graph image to selectively
hide certain edges in the final graph. This can be used to hide edges
between paris of metabolites that might have many connections in the
final graph images. Typical examples of these pairs include ATP and
ADP, NAD and NADH, etc. To use this option first a tab separated
table containing the metabolite pairs to hide must be made:

.. code-block:: shell

    atp[c]  adp[c]
    h2o[c]  h[c]
    nad[c]  nadh[c]

This file can then be used with the ``vis`` command through
the ``--hide-edges`` argument:

.. code-block:: shell

    $ psamm-model vis --hide-edges {path to edges file}


Graph Image Customization
~~~~~~~~~~~~~~~~~~~~~~~~~

By default, the reaction and metabolite nodes in a graph will only
show the reaction or metabolite IDs, but the final graphs output by the command
can be customized to include additional reaction metabolite information
that is present in the model. This additional information will be shown
directly on the reaction or metabolite nodes in the graph.
This can be done through using the ``--rxn-detail`` and ``--cpd-detail``
options. These options can be used followed by a space separated list
of properties to include. For example, the following command could be used to
show additional information of reactions and compounds:

.. code-block:: shell

    $ psamm-model vis --cpd-detail id formula charge \
      --rxn-detail id name equation


The reaction and metabolite nodes can be further customized by specifying the
filling color of nodes. This can be done by providing a two-column file that
contains reaction or metabolite IDs (with compartment) and hex color codes:

.. code-block:: shell

    ACONTa  #c4a0ef
    succ[c] #10ea88
    FUM #c4a0ef
    ....

This file can be used to color the nodes on the graph
through the ``--color`` option:

.. code-block:: shell

    $ psamm-model vis --color {path to color table}


The graph image can be simplified through the use of the ``--combine`` option.
By default, the combine level is 0. The graph generated from using
combine level 0 will have one reaction node for each reactant product
pair within a reaction. This can result in having many sets of
substrates/reaction/product nodes within the graph image, depending on how
many substrates and products are present in a metabolic reaction. Using
the combine level 1 option will condense the reaction nodes down so that
there is only one reaction node per reaction, with each reaction node having
connections to all reactants and products of that reaction. The combine level
2 option will condense the graph in a different way. With this option
the graph is condensed based on shared reactant/product pairs between
different reactions. If two separate reactions contain a common
reactant/product pair, such as ATP/ADP pair, then the nodes for those
condensed into one combined node.

.. code-block:: shell

    $ psamm-model vis --combine {0,1,2}


In some cases the exported image contains many small isolated component
that may cause the image too wide and hard to view. ``--array`` option
can be used in this cases to get a better layout. This option is
followed by an integer that is larger than 0, which indicates how many
isolated components will be placed per row. The command looks like the
following:

.. code-block:: shell

    $ psamm-model vis --array {integer that is larger than 0}

The final graph image can also be modified to show the reactions and
metabolites in different compartments based on the compartment information
provided in the model's reactions. This can be done through using the
``--compartment`` option:

.. code-block:: shell

    $ psamm-model vis --compartment

Users can specify name of output through  ``--output`` option. By default,
output will be named "reactions.dot", "reactions.nodes.tsv",
"reactions.edges.tsv":

.. code-block:: shell

    $ psamm-model vis --output Ecolicore

The output files will be named "Ecolicore.dot", "Ecolicore.nodes.tsv",
"Ecolicore.edges.tsv".


The image file produced from the ``vis`` will be automatically sized by
the `Graphviz` programs used to generate the image file. If a specific
size is desired the ``--image-size`` argument can be used to supple a
width and height in inches that the final image file should be. For
example to generate a graph that will be made into a 5" width by 10"
height image the following command can be used:

.. code-block:: shell

    $ psamm-model vis --image {format (pdf, eps, svg, etc.)} --image-size 5 10


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
specific compound. If more than one compound identifier is given
(comma-separated) this will find reactions that include all of the given
compounds.

PSAMM-SBML-Model
----------------
`PSAMM` normally takes a model in the `YAML` format as input. To deal with
models that are in the `SBML` `PSAMM` includes various programs that allow
users to convert these models to the `YAML` format. One additional option
for dealing with models in the `SBML` format is using the `psamm-sbml-model`
function. This function can be used to run any command normally accessed
through `psamm-model` but with an `SBML` model as the input. To use this
command the `SBML` model file needs to be specified first followed
by the commands:

.. code-block:: shell

    $ psamm-sbml-model {model.xml} fba

Console (``console``)
---------------------

This command will start a Python session where the model has been loaded into
the corresponding Python object representation.

.. code-block:: shell

    $ psamm-model console


Psammotate (``psammotate``)
---------------------------

Given a reciprocal best hits file, will generate a draft model based on an
template based on gene associations provided by the template file/reference
file gene mapping. Draft model will contain all relevant model components
in yaml format.

.. code-block:: shell

    $ psamm-model psammotate

To generate a draft model, a reciprocal best hits file must be specified that
maps the draft model genes to a template model using the ``--rbh`` option.
Within this file, you must specify the integer of the column that contains the
target mapping and the column that contains the template mapping (both indexed
from 1) using the ``--target`` and ``--template`` options, respectively.

.. code-block:: shell

    $ psamm-model psammotate --rbh gene_mapping.tsv --template 1 --target 2

Typically, this program retains reactions that have no gene mappings; however, if
you want to drop reactions that do not have gene associations, you must specify
the ``--ignore-na`` option.

.. code-block:: shell

    $ psamm-model psammotate --rbh gene_mapping.tsv --template 1 --target 2 --ignore-na

You can also specify an output directory for all of the yaml file output using
the ``--output`` option.

.. code-block:: shell

    $ psamm-model psammotate --rbh gene_mapping.tsv --template 1 --target 2 --output out

GIMME (``gimme``)
-----------------

This command allows you to subset a metabolic model based on gene expression
data. The expression data for filtering may be in any normalized format (TPM,
RPK, etc.), but the threshold value supplied to gimme must be appropriate for the
input data. Gimme functions through gene inactivation and will not "express" genes
that do not meet the specified expression threshold. Expression thresholds can
be specified using the ``--expression-threshold`` argument and a file that maps
genes in the model to their expression can be provided using the option
``--transcriptome-file``.

.. code-block:: shell

    $ psamm-model gimme --transcriptome-file file.tsv --expression-threshold 1

The gimme command may also specify an argument that will retain any reactions
required in order to maintain a specific biomass threshold. This threshold may
be specified using the ``--biomass-threshold`` option.

.. code-block:: shell

    $ psamm-model gimme --transcriptome-file file.tsv --expression-threshold 1 --biomass-threshold 1

You can specify a directory to output the subset model that will create all
yaml files for the new, subset, model in this directory. This location can
be specified using the ``--export-model`` argument.

.. code-block:: shell

    $ psamm-model --transcriptome-file file.tsv --expression-threshold 1 --export-model ./gimme_out/

.. _commands-tmfa:

TMFA (``tmfa``)
---------------

This command can be used to run growth simulations using the TMFA algorithm
as detailed in [Henry07]_. This method simulates fluxes in a metabolic model
while incorporating metabolite concentrations and thermodynamic constraints.
This is done through the incorporation of the gibbs free energy equation along
with additional constraints to link reactions fluxes to thermodynamic
feasibility.

This simulation method requires multiple input files to be prepared beforehand.
The following section is going to detail these input files and their formats.

TMFA Input Files
~~~~~~~~~~~~~~~~

Gibbs Free Energy Files
^^^^^^^^^^^^^^^^^^^^^^^

The ``tmfa`` method in psamm requires standard gibbs free energy of reaction
values to be supplied as an input. These values could be obtained from a
database or predicted based on metabolite structures and reaction equations.
These values must be in kilojoules per mol (kJ/mol). The input file is formatted
as a three column, tab separated table, consisting of reaction IDs, standard
gibbs free energy values, and gibbs free energy uncertainty values for the
predicted values. Any reactions that do not have associated gibbs
free energy predictions need to either be included in the lumped
reactions or in the excluded reactions.

.. code-block:: shell

    rxn1  -5.8  1.2
    rxn2  12.1  2.5
    rxn3  -0.1  0.2
    ....

Excluded Reactions File
^^^^^^^^^^^^^^^^^^^^^^^
Reactions that cannot be thermodynamically constrained need to be included in
the list of excluded reactions. These reactions typically consist of reactions
that are artificial (for example, biomass equations) or ones for which the
gibbs free energy is unknown. This input file just consists of reaction IDs
with each line containing just one ID. Note that exchange reactions will be
automatically excluded and do not need to be included in this file.

.. code-block:: shell

    biomass
    rxn6
    rxn7
    ....

Lumped Reactions File
^^^^^^^^^^^^^^^^^^^^^^^
The lumped reactions file is where lumped reactions are defined for the
TMFA problem. These reactions are summary reactions that can be used to
constrain groups of reactions, when some of the reactions in the groups
have unknown gibbs free energy values. More details on this concept can be
found in the original TMFA publication [Henry07]_. This file consists of
3 columns where the first column is a lumped reaction ID, the second column
consists of a comma separated list of reaction IDs and directions (indicated
by a 1 for forward and -1 for reverse), and a lump reaction equation. Any Lumped
reactions need to also have their gibbs free energy values assigned in the
gibbs free energy input file. The reactions that are included in the lumped
reactions do not need to be included in the excluded reactions file.

.. code-block:: shell

    Lump1  GLYBt2r:-1,GLYBabc:1    atp[c] + h2o[c] <=> adp[c] + h[e] + pi[c]
    Lump2  ADMDCr:1,SPMS:1,METAT:1 atp[c] + h2o[c] + met-L[c] + ptrc[c] <=> 5mta[c] + co2[c] + pi[c] + ppi[c] + spmd[c]


Transporter Parameters
^^^^^^^^^^^^^^^^^^^^^^^
Gibbs free energy values for transporters have to also account for the transport
of charge and pH across compartments in a model. To allow these calculations to
be made each transporter in the model needs to have an associated gibbs free
energy values defined in the gibbs free energy input file along with
the counts of protons and charge transported across the membrane defined in the
transport parameter input file. This file is formatted as three columns with
the first being reaction ID for the transporter reactions, the second being the
net charge transported from the outside compartment, and the third being the
number of protons transported from out to the in compartment for that reaction.

As an example, if we had the following reaction equation where H is a proton with
a +1 charge and X is a compound with 0 charge.

.. code-block:: yaml

    - id: transX
      equation: X[e] + H[e] <=> X[c] + H[c]

This equation would have a net charge transported of +1 and the net number of
protons transported from out to in of 1. This would mean that the transporter
parameter line for this reaction would be:

.. code-block:: shell

    transX  1 1

Another case could be if the compounds were both charged. In this case the
net charge needs to be used. For example, in the following reaction where
the compound Y is now going to have a charge of -1.

.. code-block:: yaml

    - id: transY
      equation: Y[e] + H[e] <=> Y[c] + H[c]

The protons transported from out to in would still be 1, but since compound Y
also has a charge in this case the net charge transported would be 0 (+1 for
the proton and -1 for Y). This would have an input line of the following in
the transport parameter file:

.. code-block:: shell

    transY  0  1


If reactions involve antiport, where one compound goes in and one goes out,
then this needs to be accounted for in these transporter values. The
values are always calculated in the direction of out to in. So, if there
was the export of 1 proton in the reaction equation, then the value for
that reaction would be -1 protons transported from out to in.

Lastly many reactions involve other components to their equations related
to energy costs for the reaction. For example, PTS transport reactions
may also involve the conversion of PEP to Pyruvate. These other compounds
are not considered in the transport parameter calculations, even if they
do have a charge. The only metabolites that are considered are the ones
that cross the membrane in the reaction equation.


Concentration Settings
^^^^^^^^^^^^^^^^^^^^^^^
The last input file that needs to be prepared is the concentrations file.
This file is used to set the concentrations of any metabolites that in
the model. By default, metabolites are allowed to vary in concentration
between 0.02 M and 1e-6 M. This file can be used to designate any
metabolites that will have non-default bounds in the model. These could
include ones where the measured concentrations fall outside of the
default range, or ones where the they have measured values that may
help constrain the model. For example, metabolites that are provided
in the media can have their concentrations set in this file. The
file is set up as a 3 column, tab-separated, table where the first
column is the metabolite ID (with compartment ID), the second column
is the lower bound of the concentration, and the third column is the
upper bound of the concentration. All concentrations are designated as
molar concentrations.

.. code-block:: shell

    cpd_a[e]    0.02  0.02
    cpd_b[e]    0.1 0.1
    cpd_e[e]    0.0004  0.0004


Configuration File
^^^^^^^^^^^^^^^^^^^
Because the ``tmfa`` function requires multiple input files to run,
the command was organized to use a central configuration file.
A template configuration file can be generated using the following
`PSAMM` command:

.. code-block:: shell

    $ psamm-model tmfa --config ./config.yaml util --generate-config


The configuration file contains the relative paths to the input
files that are detailed above, as well as a few other TMFA specific
parameters. Note: The relative paths in this file have to be set
based on where you are running the command from, not based on
where the configuration file is located. Additionally, any
files which are not needed for a specific model can be left out
of this configuration file.

.. code-block:: yaml

    deltaG: ./gibbs-with-uncertainty.tsv
    exclude: ./exclude.tsv
    transporters: ./transport-parameters.tsv
    concentrations: ./concentrations.tsv
    lumped_reactions: ./lumped-reactions.tsv
    proton-in: h[c]
    proton-out: h[e]
    proton-other:
     - h[p]
     - h[mito]
    water:
      - h2o[c]
      - h2o[e]
      - h2o[mito]
      - h2o[p]

In addition to the relative paths to the input files, the configuration
file is also where metabolites that are considered special cases
are designated. Specifically, water and protons are not considered as
part of the Keq portion of the gibbs free energy calculations. Becuase
of their special properties, they are considered separately from the
other metabolites. Protons are also involved in the calculation of the
gibbs free energy for transporter reactions. As such the proton ID for
the internal compartment protons and external compartment protons
need to be designated in this file as well.


TMFA Simulation Command Line Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The ``tmfa`` command has 4 general settings that can be applied to any
simulation. These settings are things that might be changed
from simulation to simulation, so they are defined through command line
parameters instead of in the configuration file.

Temperature
^^^^^^^^^^^^
The temperature used in the calculation of the gibbs free energy of the
reactions can be set through the ``--temp`` command line parameter. This
temperature is provided in Celsius.

.. code-block:: shell

    $ psamm-model tmfa --config ./config.yaml --temp 27 simulation


Biomass Flux
^^^^^^^^^^^^^
The biomass flux can be fixed at a certain value in the simulation using
the ``--threshold`` command line option. By default, the biomass is fixed
at the maximum biomass value for the model.

.. code-block:: shell

    $ psamm-model tmfa --config ./config.yaml --threshold 0.1 simulation

Reaction Reversibility
^^^^^^^^^^^^^^^^^^^^^^^
All of the reactions in a model can be allowed to be reversible, and have their
actual reversibility in the simulation only determined by the thermodynamic
constraints in the simulation. This follows a variation of TMFA as detailed in
the paper [Hamilton13]_. This can be enabled for a simulation using the
``--hamilton`` command line argument.

.. code-block:: shell

    $ psamm-model tmfa --config ./config.yaml --hamilton simulation

Error Estimates
^^^^^^^^^^^^^^^^
The gibbs free energy error estimates can be incorporated into the ``tmfa``
simulations by using the ``--err`` command line argument. This argument will
allow for the gibbs free energy values to vary slightly based on the provided
uncertainty estimates in the input file.

.. code-block:: shell

    $ psamm-model tmfa --config ./config.yaml --err simulation


Running Simulations with TMFA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Variability Analysis
^^^^^^^^^^^^^^^^^^^^^
The default simulation method with TMFA calculates the lower and upper bounds
for each variable in the TMFA problem. this produces a four-column output that
gives the variable type, variable ID, lower bound, and upper bound.
This variability analysis is done for the metabolite concentrations, reaction
fluxes, reaction gibbs free energy values, and binary constraint variables for
the reactions. This analysis can be run using the following command:

.. code-block:: shell

    $ psamm-model tmfa --config ./config.yaml simulation


This will produce an output like this:

.. code-blocks:: shell

    Flux	CS	0.18355092011356122	0.18355092011362384
    DGR	CS	-42.528916818320724	-1.0000000010279564e-06
    Zi	CS	1.0	1.0
    CONC	icit_c[c]	9.999e-06	0.007
    CONC	co2_c[c]	9.999e-06	9.999e-05
    CONC	2pg_c[c]	9.999e-06	0.0127

Single Solution Simulations
^^^^^^^^^^^^^^^^^^^^^^^^^^
The TMFA simulations can also be run to produce single, arbitrary solutions
instead of doing the variability analysis. This can be run to produce a
single solution using just the TMFA constraints ``fba``, by doing an additional
L1 minimization of the fluxes ``l1min``, or to produce a random solution from
the solution space ``random``. (NOTE: the random solution is experimental and
is not guaranteed to be completely random).

.. code-block:: shell

    $ psamm-model tmfa --config ./config.yaml simulation --single-solution fba

    $ psamm-model tmfa --config ./config.yaml simulation --single-solution l1min

    $ psamm-model tmfa --config ./config.yaml simulation --single-solution random

Randomsparse with TMFA Constraints
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Random reaction or gene deletions can be performed while accounting for the
thermodynamic constraints through using the ``--randomsparse`` or
``--randomsparse_genes`` command line arguments. These functions will
perform random deletions of reactions or genes in the model in the same way
as the PSAMM ``randomsparse`` command, but while accounting for the
thermodynamic constraints in the model.

.. code-block:: shell

    $ psamm-model tmfa --config ./config.yaml simulation --randomsparse

    $ psamm-model tmfa --config ./config.yaml simulation --randomsparse_genes


Other TMFA Utility Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The ``tmfa`` utility functions include the ``--generate-config`` command
, which was detailed above, and a testing command ``--random-addtion`` which
can be run to test which thermodynamic constraints might be causing
a model to fail to produce biomass. In some cases when the simulations settings
are being set up or when the parameters are changed (for example, different temps
or concentrations) the ``tmfa`` simulation might not produce biomass.
Unfortunately, it can be difficult to isolate which constraints might cause
this problem without extensive manual investigation. To help with this process
the ``--random-addition`` utility function was developed. This function will
randomly add reaction thermodynamic constraints to the model and test which
constraints cause the model to fall below the set biomass threshold. This can
be useful for identifying potentially problematic or over constrained reactions
in the model.

.. code-block:: shell

    $ psamm-model tmfa --config ./config.yaml util --random-addition
