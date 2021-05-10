v1.1.2 (2021-05-10)  
------------------
  
- Adds the --phin and --phout options to the `tmfa` command
- Adds descriptive counts of the number of reactions and compounds
  curated in the _manual curation_ option of `modelmapping`


v1.1.1 (2020-12-18)
------------------

- Fix bug where the `psamm modelmapping translate_id` command requires 
  the undefined parameter `--dest-model`


v1.1 (2020-12-11)
------------------

- Adds the new `vis` command for visualizing metabolic pathways and 
  networks
- Adds the new `tmfa` command to implement the approach as described 
  in Henry et al., 2007 (PMID: 17172310) and Hamilton et al., 2013 
  (PMID: 23870272).
- Performance updates to the `model-mapping` command
- Updates to export full models from `gimme` and `psammotate`
- Added the --fva option to the `robustness` command


v1.0 (2020-03-03)
------------------

- Adds the new `modelmapping` command for mapping reactions and
  compounds ids between different GEMs.
- Adds the new `psammotate` command for generating draft GEMs
  using a template model and a mapping of orthologous genes to
  another organism.
- Adds the new `gimme` command which implements the GIMME algorithm
  as described in Becker and Palsson, 2008 (PMID: 18483554).
- Updates the `search` command to allow for searching of strings
  within any reaction or compound property.
- Updates to the `tableexport` and `excelexport` commands to allow
  for the export of additional gene and reaction information.
- Adds new section to the tutorial to detail how to use the
  `findprimarypairs` commands.
- Renamed `duplicatescheck` command to `dupcheck`.

v0.31 (2017-06-30)
------------------

- The `psamm-import` tool has been moved from the `psamm-import` package to
  the main PSAMM package. This means that to import SBML files the
  `psamm-import` package is no longer needed. To use the model-specific Excel
  importers, the `psamm-import` package is still needed. With this release
  of PSAMM, the `psamm-import` package should be updated to at least 0.16.
- The tutorial was updated with additional sections on using gap-filling
  procedures on models.

v0.30 (2017-06-23)
------------------

- Adds the new command `primarypairs` for predicting reactant/product element
  transfers using the new FindPrimaryPairs method as well as the MapMaker
  method.
- A new option has been added to the `genedelete` command which allows use of
  _minimization of metabolic adjustments_ to simulate biomass production for
  gene knockouts.
- Fixes a bug where the epsilon parameter was accidentally ignored by
  `fastgapfill`.
- Fixes a bug where the `psamm-sbml-model` command did not ignore boundary
  species. With this change, the boundary species are also ignored by default
  when using the API to read SBML models.
- Fixes a performance issues with gap-filling that made the `gapfill` and
  `fastgapfill` commands take much longer time to run on large models than
  necessary.

v0.29 (2017-04-18)
------------------

- [The tutorial](https://psamm.readthedocs.io/en/stable/tutorial.html) in the
  PSAMM documentation has been updated and expanded to include additional
  information on using PSAMM for model curation and constraint-based analyses.
- The experimental command `psamm-sbml-model` was added which makes it possible
  to run any command from `psamm-model` (e.g. `fba`, `robustness`, etc.)
  directly on an SBML file. For now this only supports SBML level 3 files with
  FBC. This provides a quick way of running basic analyses on SBML files. We
  still recommend importing the SBML file to YAML format with `psamm-import`
  for anyone wishing to make changes to a model.
- Fixes access to charge parameter parsed from SBML files. The charge is now
  correctly imported with `psamm-import`.
- Fixes import of compartments from SBML files. The empty boundary compartments
  are now no longer included in the import.
- Fixes bug in writing the reaction flux limits sheet of the `excelexport`
  command.
- The `console` command was changed to only provide the `model` variable since
  the metabolic model can easily be created.

v0.28 (2017-03-03)
------------------

- The YAML model format now allows users to specify compartment information and
  compartment boundaries in the `model.yaml` file. See the file format
  documentation for more information.
- The `media` key in the `model.yaml` has changed name to `exchange` to
  reflect the fact that not only uptake exchange must be defined here. The
  `media` key is still supported but has been deprecated.
- The gap-filling command `gapfill` and `fastgapfill` now use the compartment
  information to determine which artificial transport and exchange reactions
  to add. This means that a model *must* specify compartments and compartment
  boundaries when using gap-filling commands.
- The `gapcheck` command now has two new methods for detecting blocked
  compounds. The new `prodcheck` is a more robust version of the GapFind check
  which was previously used. The new `sinkcheck` method will find compounds
  that cannot be produced in excess. This can find some additional blocked
  compounds that were not detected by the other methods.
- The `gapcheck` command now reports blocked compounds in the extracellular
  space. Previously, these compounds were excluded. An option is available to
  switch back the old behavior of excluding these from the final output.
- The `gapcheck` command now has an option to run the check with unrestricted
  exchange reactions.
- The `gapfill` command can now be run without implicit sinks. This makes it
  possible to use this command to solve additional model gaps. It is still
  recommended to first solve gaps using implicit sinks, then later disable
  implicit sinks when all other gaps have been closed.
- The `gapfill` command now has an option to enable the bounds expansion
  proposals (e.g. making irreversible reactions reversible). By default this
  option is now off.
- The `fastgapfill` has improved output that contains less superfluous
  information. The output format is now identical to the `gapfill` command.
  The `fastgapfill` also no longer runs an FBA on the induced model since this
  caused some confusion.
- Added new command `checkduplicates` which detects whether the model has
  multiple reactions with the same (or similar) reaction equation.
- The `sbmlexport` command now allows the user to specify a file path. The
  command can also optionally output the SBML file in a more readable format
  with an option.
- Fixed support for the latest CPLEX release 12.7. A change in their API made
  PSAMM incompatible with the 12.7 release. This is now fixed.
- We now officially support Python 3.5 and Python 3.6.

v0.27 (2016-12-23)
------------------

- When exporting a model with `excelexport`, `sbmlexport` or `tableexport`, the
  current Git commit ID of the model is now part of the export.
- The `gapfill` command has been split into two separate commands. The new
  command `gapcheck` only reports which compounds are blocked without
  trying to fill gaps. The `gapfill` command now only reports the suggested
  reactions for gap-filling.
- The `tableexport` command has been made more robust when the model properties
  contain special characters or structured data. Strings containing tabs or
  newline characters are now quoted in the output.
- Improved error messages for multiple commands when reactions with forced
  flux (e.g. ATP maintenance reactions) impose impossible constraints.

v0.26 (2016-11-17)
------------------

- All commands that perform loop removal now use the option `--loop-removal` to
  set the type of loop removal. The `--tfba` option is no longer available for
  the `fva`, `fluxcheck` and `randomsparse` commands, instead
  `--loop-removal=tfba` should be used.
- All simulation commands that perform a maximization of an objective now take
  an `--objective` option that can be used to override the `biomass` reaction
  set in the model. The commands `fba` and `fva` no longer allow positional
  arguments to set the objective, instead the `--objective` option should be
  used.
- All user data in the model is now exported with the `sbmlexport` command.
  User data that cannot be translated into a standard machine-readable SBML
  form will instead be stored in the SBML notes section.
- The output from the `gapfill` command now shows the complete induced model
  (like `fastgapfill`) and also lists the compound IDs of the blocked
  compounds.
- The `gapfill` command will now propose the removal of flux bounds if such a
  modification can unblock the model.
- The `gapfill` command now excludes the biomass reaction from being modified
  in the induced model. Previously the biomass reaction was sometimes reversed
  in the induced model but this is usually not a desired solution.
- The `gapfill` command now allows the user to specify penalties on added or
  modified reactions (like `fastgapfill`).
- The `gapfill` command was changed so the `--epsilon` option now
  specifies the threshold for non-zero reaction fluxes (like the `--epsilon`
  in other PSAMM commands).
- The `gapfill` command was changed to be more robust when handling small
  epsilon values. It will now lower the solver thresholds if necessary and will
  generate a warning when epsilon is too low. This fixes an issue where
  `gapfill` would previously fail or generate incorrect results with some
  models.
- Fixed: The `--exclude` option in the `formulacheck` command did not work
  correctly.

v0.25 (2016-10-28)
------------------

- Fixed an error parsing decimal values in reactions which resulted in a
  failure to run FBA and other analyses on certain models.
- The `fastgapfill` command no longer tries to unblock exchange reactions by
  default. Only internal reactions will be unblocked by default.
- The `fastgapfill` command has a new `--subset` option to explicitly specify
  set of reactions to unblock. This means that the command can now be used to
  unblock a specific reaction.
- The weight options on `fastgapfill` have changed name to `--db-penalty`,
  `--tp-penalty` and `--ex-penalty` for consistency with the existing
  `--penalty` option.
- Fixed an error in `gapfill` that in some cases would result in a compound
  incorrectly marked as non-blocked.
- The `sbmlexport` command now follows the FBCv2 specification for writing
  flux bounds, biomass reaction, gene products and various other properties to
  SBML files.
- The `sbmlexport` command now uses the same ID for compartment IDs as used
  in the YAML files.
- The order of compounds in the output from commands now reflects the order
  of compounds in the reactions as specified in the model files.
- Experimental support for solving MILP problems with GLPK has been activated.

v0.24 (2016-08-23)
------------------
- When specifying flux bounds in media and limits, the `fixed` key can now be
  used when lower and upper limits are the same.
- New column in output of `excelexport` and `tableexport` commands to
  indicate if reactions and compounds are in the model.
- Zero mass compounds are now omitted from the output of the `massconsistency`
  command.
- When exporting SBML file using the `sbmlexport` command, the exported
  compounds and reactions now have IDs that are based on the YAML model IDs.
  Characters that are not allowed in SBML IDs are transformed in a way that is
  compatible with COBRA.

v0.23 (2016-07-26)
------------------

- Fix a bug where no output of the `randomsparse` command was produced.
- Make Cplex interface in PSAMM compatible with earlier versions
  (12.6 and 12.6.1) again.

v0.22 (2016-07-01)
------------------

- Better unicode handling in commands.
- When running the `gapfill` command the epsilon parameter can now be
  specified on the command line.
- When parsing reaction and compound entities from the YAML files, produce
  better error messages when IDs are invalid.
- Work around a bug in Cplex that in rare causes a segmentation fault when a
  linear programming problem is solved repeatedly.
- API: Add `fastgapfill` module which allows access to run the fastGapFill
  algorithm.
- API: Add `randomsparse` module which allows access to generate a random
  minimal model which satisfies the flux threshold of the objective reaction.

v0.21 (2016-06-09)
------------------

- Add `genedelete` command to allow users to delete one or more genes and
  perform a viability check on the model after all related reactions are
  deleted.
- Add `balancecheck` module which allows API access to charge balance and
  formula balance checks.
- When a compound in the extracellular space doesn't have an exchange reaction,
  a warning would be provided so that the user may add the compound to the
  medium.
- If a compartment is not given for a medium, it will now be assumed to have
  the extracellular compartment.
- When using Gurobi and Cplex, the default optimality tolerance and feasibility
  tolerance has been decreased to 1e-9.
- Fixed a bug where the reaction IDs are not printed properly in the result of
  `chargecheck` command.

v0.20 (2016-03-10)
------------------

- Added experimental support for GLPK solver. MILP problems are not yet
  supported with this solver. GLPK also appears to have some issues with
  `fastgapfill`.
- The `gapfill` command can now take a list of compounds on the command line
  that it will try to unblock. If a list of compounds is given, the command
  will not run GapFind but instead only use those compounds.
- Remove the assumption that the extracellular compartment is always called
  `e`. The user can now specify the name of the extracellular compartment with
  the option `extracellular` in `model.yaml`.
- In a previous release, the code interfacing with Cplex was updated and is now
  using an interface in Cplex that was introduced in version 12.6.2. The
  documentation now makes it clear that at least version 12.6.2 is required.
- Update YAML format documentation on the model definition table format.
- Work around issue with `pkg_resources` that resulted in import errors when
  running from IPython.

v0.19 (2016-02-01)
------------------

- When using Gurobi, the option `--solver threads=X` can now be used to specify
  the maximum number of threads that Gurobi can use.
- The log messages from external libraries are more clearly marked as such.
  In particular, there should now be less confusion about the origin of the
  log messages from Cplex.
- Internally, the LP solvers now support quadratic objectives. This will be
  used for various commands in the future.
- Fix an error where an empty reaction would internally be detected as an
  exchange reaction.
- Fix a bug where the compounds in the extracellular compartment were not
  correctly detected by `gapfill`.
- Update documentation with information on how to cite the PSAMM publication.

v0.18 (2016-01-04)
------------------

- Several commands now support parallelization with the `--parallel` option
  (`fva`, `fluxcheck`, `fluxcoupling`, `robustness`).
- A more robust reaction parser is now used to parse reaction equations in
  YAML files. This also means that quoting compound names with pipes (`|`) is
  now optional.

v0.17 (2015-12-07)
------------------

- When loading native models, PSAMM now uses the PyYAML safe loader and also
  uses the optimized CSafeLoader if present. This speeds up the start time of
  commands.
- Various additional optimizations to model loading have been added. This
  speeds up the start time of some commands.
- The `fba` command now shows the genes associated with each reaction for a
  quick overview of which genes influence the flux solution.
- The `sbmlexport` command now properly exports gene association information.
- All commands better handle output that contains unicode characters. In
  previous versions this would often fail when using Python 2.

v0.16 (2015-12-01)
------------------

- Add an option to `randomsparse` to perform the deletion based on genes
  instead of reactions. This uses the gene association expression defined in
  the reaction property `genes`.
- Add threshold option to `fva` command.
- Fix bugs in `gapfill` that resulted in the procedure not detecting reactions
  that could be reversed, and sometimes failing to find a result at all.
- Add epsilon option to `chargecheck` and ignore charge imbalances below the
  epsilon value.
- Allow the `search` command to find reactions containing a specific compound
  even when the compartment is not specified.
- Output more information is the result of the `search` command.
- Improved handling of flux bounds at infinity (e.g. with
  `default_flux_limit = .inf` in `model.yaml`).

v0.15 (2015-10-16)
------------------

- Add support for reading flux bounds and objectives from SBML files that are
  using the FBC extension.
- Add a tutorial to the documentation at <https://psamm.rtfd.org/>.
- Add command `tableexport` to export various parts of the model as a TSV file.
- Add command `excelexport` to export all parts of the model as an Excel file.
- Allow various parameters that take a reaction as an argument to also be able
  to take a list of reactions from a file, using the `@` prefix. For example,
  given a file `r.txt` with each reaction ID on a separate line, the reactions
  can be excluded from the `masscheck` command by specifying
  `--exclude @r.txt`.
- Allow reactions to be excluded from the `formulacheck` and `chargecheck`
  commands.

v0.14 (2015-10-05)
------------------

- Split the `masscheck` command into two parts. The compound check is run be
  default or when `--type compound` is specified. The reaction check is run
  when `--type reaction` is specified. This also changes the output of the
  compound list to be TSV formatted.
- Change the default method of finding flux inconsistent reactions to simply
  use FVA without constraints to determine whether reactions can take a
  non-zero flux. This requires more LP optimizations to run but it turns out
  to be faster in practice. To enable the old behavior where the number of LP
  problems to solve is reduced, use `--reduce-lp`.
- Disable tFBA by default in the FBA performed as part of running
  `fastgapfill`.
- Return non-zero from the `psamm-list-lpsolvers` when no solver is available.
- Report time to solve most commands excluding the time it takes to load the
  model.
- Improve stability when using thermodynamic constraints. This means that
  commands using thermodynamic constraints that previously failed with some
  models will now work.
- Speed up changing the objective when using Cplex. This significantly speeds
  up commands that reuse LP problem instances with different objectives (e.g.
  `fva` and `fluxcheck`).
- Speed up fastcore algorithms (i.e. `fluxcheck --fastcore` and `fastgapfill`)
  by reusing the LP problem instances.
- Propagate user aborts from Cplex to Python by raising `KeyboardInterrupt`
  when a user abort is detected. This fixes a problem where a user abort would
  result in a `FluxBalanceError`.
- Improve unit tests of commands, the `native` datasource module,
  `fluxanalysis`, and various other parts of the software.

v0.13 (2015-09-16)
------------------

- Change parsing of medium definitions to combine all parts of the medium into
  one final medium used for simulations. Previously, the entries in the
  `media` key would be considered as separate media and only the first entry
  would be used.
- Fix parsing of limits table format to make it consistent with the medium
  table format.
- Change the name of option `--reaction` to `--objective` for commands that
  take an optional reaction flux to maximize.
- Change the default of commands that use thermodynamic constraints on FBA
  problems to _not_ apply the additional constraints by default. The commands
  have an option to enable thermodynamic constraints.
- Change the `--no-tfba` option of the commands `fba` to `robustness` to
  `--loop-removal`. The loop removal option allows the L1 minimization method
  to be selected as a loop removal method in addition to the thermodynamic
  constraints.
- Add option to `robustness` command to print all fluxes. By default only the
  objective is shown.
- In `robustness`, obtain the lower and upper bounds of the varying reaction
  by minimizing and maximizing the reaction flux. This avoids trying to solve
  many LP problems that are infeasible and instead provides more samples from
  the flux range that is known to be feasible.
- Fix a bug in `robustness` that caused thermodynamic constraints to never be
  enabled properly.
- Fix bug where parsing medium from table would fail.
- Produce a warning when a compound is referenced from a reaction but not
  defined in the list of compounds.
- During SBML import produce a warning when `_b`-suffix compounds are
  converted to boundary compounds in non-strict mode.
- Improve the efficiency of creating long linear expressions for the LP
  solver interface. This should speed up certain procedures that use long
  expression as objectives or constraints.
- Improve the efficiency of resolving variable names in Cplex LP problems by
  using column indices directly. This significantly speeds up procedures using
  many variables when using Cplex as the solver.

v0.12.1 (2015-09-01)
--------------------

- Fix tiny bug in setup.py resulting in failure to upload the 0.12 package to
  PyPI.

v0.12 (2015-09-01)
------------------

- Add support for the Gurobi LP solver (Python 2.7 only).
- Fix a bug in the `fastgapfill` command that caused the biomass reaction to
  never be found in the model.
- Add command line tool `psamm-list-lpsolvers` which will list the available
  solvers, and report their attributes and whether they can be successfully
  loaded.
- Change the way the threshold parameter on the `randomsparse` command is
  parsed. A threshold relative to the maximum flux is now given using
  percentage notation (e.g. `95%`) while an absolute flux threshold is given as
  a number (e.g `1.34`).
- Log the model name at the beginning of every command. If the model exists in
  a Git repository, the current Git commit ID is also reported.
- Produce warnings when encountering invalid SBML constructs using non-strict
  parsing. Previously, some of these could only be noticed when using the
  strict mode.
- Improve error reports from commands. This means that usage information is now
  printed when a command fails to parse an argument.
- Improve the API in the `fluxanalysis` module to allow commands to more
  easily reuse LP problems. This improves the speed of robustness analysis in
  particular.
- Improve and expand the install instructions in the documentation.
- Reorganize implementations of commands into the `psamm.commands` package. The
  commands are now discovered through the entry point mechanism from
  `setuptools` allowing any package to install additional commands into the
  PSAMM user interface.
- Use a simpler output format for log messages unless the environment variable
  `PSAMM_DEBUG` is set. If set, this variable determines the logging level
  (e.g. `PSAMM_DEBUG=debug` enables more messages).
- Include Python 3.3 in the tox test suite.

v0.11 (2015-08-06)
------------------

- Add new flux coupling analysis command and API module.
- Add parser for KEGG reaction files.
- Export COBRA-compatible flux constraints when exporting to SBML.
- Improve the command descriptions in the command line interface.
- Extend the search command to report information on which file the
  compound/reaction was parsed from.
- Add support for Python 3.4.
- Fix Fastcore unit tests that failed using an exact solver.
- Change unit tests to allow covering solvers other than Cplex.
- Optionally use tox for unit test management. This allows independently
  testing with multiple versions of Python and with different solvers.
- Adapt code to conform to PEP-8.

v0.10.2 (2015-05-22)
--------------------

- Various minor documentation and packaging updates. PSAMM is now available
  through [PyPI](https://pypi.python.org/pypi/psamm) and the documentation has
  been updated to reflect this. The documentation is available through
  [Read the Docs](https://psamm.readthedocs.org/). Lastly, the test suite is
  automatically run by [Travis CI](https://travis-ci.org/zhanglab/psamm).

v0.10.1 (2015-05-21)
--------------------

- Update README with new repository names and improved install instructions.
- docs: Add improved install instructions based on README.

v0.10 (2015-05-14)
------------------

- This software is now GPLv3 licensed. A copy of the license is included in
  [LICENSE](LICENSE).
- Allow setting the default flux limit in `model.yaml`. Previously a limit
  of 1000 units was used, and this value is still used if not specified.
- Allow setting the reaction name of the implicit exchange reactions
  specified in the medium definition.
- sbml: Change the SBML writer to avoid negative values for reaction
  stoichiometry. Some software packages do not handle negative values
  correctly when loading SBML files.
- command: Add a new option to `fluxcheck` where the restrictions imposed on
  the exchange reactions are removed before the consistency check.
- cplex: Use numerical emphasis mode by default.

v0.9 (2015-05-07)
-----------------

- Add methods the internal metabolic model representation to provide the
  compartments. This is used by commands instead of hardcoding specific
  model compartments when running `gapfill` or `fastgapfill`.
- docs: Update documentation on `psamm-model` commands.
- docs: Update information on installing the linear programming solvers.

v0.8 (2015-04-30)
-----------------

- The name of the project (in `setup.py`) changed name to `psamm`.
- The name of the main package changed from `metnet` to `psamm`.
- Remove unused scripts related to obsolete GAMS modeling.
- Assume TSV format for `.tsv` compound files. This removes the need for
  explicitly specifying `format: tsv` in the include.
- By default FVA and the flux consistency check will apply thermodynamic
  constraints. A command line option was added to go back to the previous
  behavior.
- Properly report compound names when multiple files are included.
- Add possibility of supplying additional solver parameters through the command
  line. Currently, the Cplex solver supports the parameters `threads` (the max
  number of threads allowed globally) and `feasibility_tolerance` (how much the
  basic variables of a model are allowed to violate their bounds).
- command: The command `randomsparse` now defaults to normal FBA without
  thermodynamic constraints. This is much faster and the additional constraints
  are guaranteed to not change the result in this case.
- docs: Enabled napoleon Sphinx extension for Google-style docstring support.
  This makes docstrings more readable while at the same time improving the
  generated docs. The `fluxanalysis` module was updated with additional
  documentation in this format.
- fluxanalysis: Change the API so that the `tfba` parameter can be given to
  functions that support thermodynamic constraints to select this mode. The
  support for thermodynamic constraints was extended to the flux consistency
  check and FVA.
- fluxanalysis: Slightly improve FVA by avoiding copying the model.
- sbml: Provide access to species charge.
- qsoptex: Fix error when calling the `status()` method of a result.
- command: Add option to see current version.

v0.7 (2015-04-23)
-----------------

- Change name of `model` script to `psamm-model`.
- native: Add YAML format for flux limits. The documentation has been updated
  to include more information on this format. This changes the `limits` key in
  `model.yaml` to a list of dicts.
- Add `--exchange` option to the `randomsparse` command to find random
  minimal sets of exchange reactions.
- Change name of `fluxconsistency` command to `fluxcheck`.
- Allow compounds to be marked as zero-mass (e.g. photons) by setting
  `zeromass: yes`. These compounds will be exempt from the mass requirements
  in the `masscheck` command.
- sbml: Provide access to model ID and name.
- sbml: Add option to skip boundary condition species when parsing.
- massconsistency: Fix bugs occurring when zero-mass compounds are specified.
- command: Log number of consistent reactions in `masscheck`.
- sbml: Fix a number of minor bugs.
- command: Fix search command when no alternative compound names are present

v0.6.1 (2015-04-20)
-------------------

- sbml: Fix bug where boundary conditions were parsed incorrectly.

v0.6 (2015-04-17)
-----------------

- Apply changes to the SBML parser in order for it to interoperate with
  `model-import`. This makes it easier to implement the SBML importer in
  `model-import`.
- Add non-strict mode to the SBML parser. This makes it possible to load
  almost-compliant SBML documents that are accepted by COBRA.
- `masscheck` command: Allow reactions to be marked as checked.
- cplex: Consider status `optimal_tolerance` to be successful.
- docs: Expand documentation on the `masscheck` command.
- docs: Change order of API documentation to `bysource`.

v0.5 (2015-04-09)
-----------------

- Add `sbmlexport` command to export current model as an SBML file.
- Add a generic interface to the linear programming solvers that delegates to
  an actual solver that is installed and has the required features. This adds
  a `--solver` option to a number of command which can be used to influence
  which solver is selected.
- Add `--epsilon` option to a number of commands that previously had the
  epsilon value hardcoded.
- Refactor functions in `fastcore` for easier use.
- docs: Extend docstring documentation of various modules.
- docs: Add DOI links for references.

v0.4 (2015-04-06)
-----------------

- Add documentation generated by Sphinx. The main contents of the
  [README](README.md) file has been moved to the new documentation.
- Generate the entry-point script using `setup.py`. This ensures that the
  package is correctly installed before the main script can be called. This
  also changes the name of the entry-point from `model.py` to `model`.
- Refactor functions in `massconsistency` for easier use.
- Add `__version__` attribute to main module.
- docs: Move references to separate section.
- docs: Fix file format documentation for medium file.
- Unit tests: Skip tests requiring a linear programming solver if Cplex is
  present.

v0.3 (2015-03-27)
-----------------

- Require reaction files to be explicitly listed in `model.yaml`.
- Add support for TSV reaction file format.
- Change format of YAML reactions (see [README](README.md) for details).
- Add tables of recognized compounds and reaction properties to
  [README](README.md).
- `masscheck` command: Automatically exclude biomass reaction from check.

v0.2 (2015-03-20)
-----------------

- Allow compounds to be specified using YAML, TSV or ModelSEED format. This
  changes the format of the `compounds` key in `model.yaml` (see
  [README](README.md) for more information).
- Allow specifying biomass reaction in `model.yaml` using the `biomass` key.
  The biomass reaction will be used by default for FBA, FVA, etc.
- Allow explicit definition of media. This can be defined using a table format
  or YAML format. See [README](README.md) for more information.
- `chargecheck`/`formulacheck` commands: Only check reactions where all
  compounds have charge/formula specified. The number of skipped reactions is
  reported separately.
- `chargecheck` command: Use charge information from model definition instead
  of requiring a separate charge table file.

v0.1 (2015-03-18)
-----------------

- Initial release.
