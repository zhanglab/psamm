# This file is part of PSAMM.
#
# PSAMM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PSAMM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PSAMM.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2014-2015  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2015  Keith Dufault-Thompson <keitht547@my.uri.edu>

"""Command line interface

Each command in the command line interface is implemented as a subclass of
:class:`Command`. Commands are automatically discovered based on this
inheritance so new commands are added by subclassing :class:`Command`.

The :func:`.main` function is the entry point of command line interface.
"""

from __future__ import print_function

import sys
import os
import argparse
import operator
import re
import logging
import random
import math
import abc

from . import __version__ as package_version
from .formula import Formula, Radical
from .gapfill import gapfind, gapfill
from .database import DictDatabase
from .metabolicmodel import MetabolicModel
from .reaction import Compound
from .datasource.native import NativeModel
from .datasource import sbml
from . import fluxanalysis, massconsistency, fastcore
from .lpsolver import generic

# Module-level logging
logger = logging.getLogger(__name__)


class Command(object):
    """Represents a command in the interface, operating on a model

    Subclasses must define name and title as class attributes. The constructor
    will be given the NativeModel and the command line namespace. The subclass
    must implement :meth:`run` to handle command execution.

    In addition, :meth:`init_parser` can be implemented as a classmethod which
    will allow the command to initialize an instance of
    :class:`argparse.ArgumentParser` as desired. The resulting argument
    namespace will be passed to the constructor.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, model, args):
        self._model = model
        self._args = args

        # Create metabolic model
        database = DictDatabase()
        for reaction in model.parse_reactions():
            if reaction.equation is not None:
                database.set_reaction(reaction.id, reaction.equation)

        media = list(model.parse_media())
        if len(media) > 1:
            logger.warning('Only the first medium will be used')
        medium = media[0] if len(media) > 0 else None

        self._mm = MetabolicModel.load_model(
            database, model.parse_model(), medium, model.parse_limits(),
            v_max=model.get_default_flux_limit())

    @classmethod
    def init_parser(cls, parser):
        """Initialize command line parser (:class:`argparse.ArgumentParser`)"""

    @abc.abstractmethod
    def run(self):
        """Execute command"""


class SolverCommandMixin(object):
    """Mixin for commands that use an LP solver

    This adds a ``--solver`` parameter to the command that the user can use to
    select a specific solver. It also adds the method :meth:`_get_solver` which
    will return a solver with the specified default requirements. The user
    requirements will override the default requirements.
    """

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--solver', action='append', type=str,
            help='Specify solver requirements (e.g. "rational=yes")')

    def __init__(self, *args, **kwargs):
        super(SolverCommandMixin, self).__init__(*args, **kwargs)
        self._solver_args = {}
        if self._args.solver is not None:
            for s in self._args.solver:
                key, value = generic.parse_solver_setting(s)
                self._solver_args[key] = value

    def _get_solver(self, **kwargs):
        """Return a new :class:`psamm.lpsolver.lp.Solver` instance"""
        solver_args = dict(kwargs)
        solver_args.update(self._solver_args)
        return generic.Solver(**solver_args)


class ChargeBalanceCommand(Command):
    """Check whether compound charge in a given database or model is balanced

    Balanced reactions are those reactions where the total charge
    is consistent on the left and right side of the reaction equation.
    Reactions that are not balanced will be printed out.
    """

    name = 'chargecheck'
    title = 'Check charge balance on a model or database'

    def run(self):
        """Run charge balance command"""

        # Load compound information
        compound_name = {}
        compound_charge = {}
        for compound in self._model.parse_compounds():
            compound_name[compound.id] = (
                compound.name if compound.name is not None else compound.id)
            if hasattr(compound, 'charge') and compound.charge is not None:
                compound_charge[compound.id] = compound.charge

        # Create a set of known charge-inconsistent reactions
        exchange = set()
        for reaction_id in self._mm.reactions:
            if self._mm.is_exchange(reaction_id):
                exchange.add(reaction_id)

        def reaction_charges(reaction_id):
            for compound, value in self._mm.get_reaction_values(reaction_id):
                charge = compound_charge.get(compound.name, float('nan'))
                yield charge * float(value)

        count = 0
        unbalanced = 0
        unchecked = 0
        for reaction in sorted(self._mm.reactions):
            if reaction in exchange:
                continue

            count += 1
            charge = sum(reaction_charges(reaction))
            if math.isnan(charge):
                logger.debug('Not checking reaction {};'
                             ' missing charge'.format(reaction))
                unchecked += 1
            elif charge != 0:
                unbalanced += 1
                rx = self._mm.get_reaction(reaction)
                rxt = rx.translated_compounds(
                    lambda x: compound_name.get(x, x))
                print('{}\t{}\t{}'.format(reaction, charge, rxt))

        logger.info('Unbalanced reactions: {}/{}'.format(unbalanced, count))
        logger.info('Unchecked reactions due to missing charge: {}/{}'.format(
            unchecked, count))


class ConsoleCommand(Command):
    """Start an interactive Python console with the given model loaded"""

    name = 'console'
    title = 'Start Python console with metabolic model loaded'

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--type', choices=('python', 'ipython', 'ipython-kernel'),
            default='python', help='type of console to open')

    def open_python(self, message, namespace):
        """Open interactive python console"""

        # Importing readline will in some cases print weird escape
        # characters to stdout. To avoid this we only import readline
        # and related packages at this point when we are certain
        # they are needed.
        from code import InteractiveConsole
        import readline, rlcompleter

        readline.set_completer(rlcompleter.Completer(namespace).complete)
        readline.parse_and_bind('tab: complete')
        console = InteractiveConsole(namespace)
        console.interact(message)

    def open_ipython(self, message, namespace):
        from IPython.terminal.embed import InteractiveShellEmbed
        console = InteractiveShellEmbed(user_ns=namespace, banner2=message)
        console()

    def open_ipython_kernel(self, message, namespace):
        from IPython import embed_kernel
        embed_kernel(local_ns=namespace)

    def run(self):
        message = ('Native model has been loaded into: "model"\n' +
                   'Metabolic model has been loaded into: "mm"')
        namespace = {'model': self._model, 'mm': self._mm}
        console_type = self._args.type

        if console_type == 'python':
            self.open_python(message, namespace)
        elif console_type == 'ipython':
            self.open_ipython(message, namespace)
        elif console_type == 'ipython-kernel':
            self.open_ipython_kernel(message, namespace)


class FastGapFillCommand(SolverCommandMixin, Command):
    """Run FastGapFill algorithm on a metabolic model"""

    name = 'fastgapfill'
    title = 'Run FastGapFill on a metabolic model'

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--penalty', metavar='file', type=argparse.FileType('r'),
            help='List of penalty scores for database reactions')
        parser.add_argument(
            '--db-weight', metavar='weight', type=float,
            help='Default weight for database reactions')
        parser.add_argument(
            '--tp-weight', metavar='weight', type=float,
            help='Default weight for transport reactions')
        parser.add_argument(
            '--ex-weight', metavar='weight', type=float,
            help='Default weight for exchange reactions')
        parser.add_argument(
            '--epsilon', type=float, help='Threshold for Fastcore',
            default=1e-5)
        parser.add_argument(
            '--no-tfba', help='Disable thermodynamic constraints on FBA',
            action='store_true')
        parser.add_argument(
            'reaction', help='Reaction to maximize', nargs='?')
        super(FastGapFillCommand, cls).init_parser(parser)

    def run(self):
        """Run FastGapFill command"""

        # Create solver
        enable_tfba = not self._args.no_tfba
        if enable_tfba:
            solver = self._get_solver(integer=True)
        else:
            solver = self._get_solver()

        # Load compound information
        compound_name = {}
        for compound in self._model.parse_compounds():
            compound_name[compound.id] = (
                compound.name if compound.name is not None else compound.id)

        epsilon = self._args.epsilon
        model_compartments = set(self._mm.compartments)

        # Add exchange and transport reactions to database
        model_complete = self._mm.copy()
        logger.info('Adding database, exchange and transport reactions')
        db_added = model_complete.add_all_database_reactions(model_compartments)
        ex_added = model_complete.add_all_exchange_reactions()
        tp_added = model_complete.add_all_transport_reactions()

        # Add penalty weights on reactions
        weights = {}
        if self._args.db_weight is not None:
            weights.update((rxnid, self._args.db_weight) for rxnid in db_added)
        if self._args.tp_weight is not None:
            weights.update((rxnid, self._args.tp_weight) for rxnid in tp_added)
        if self._args.ex_weight is not None:
            weights.update((rxnid, self._args.ex_weight) for rxnid in ex_added)

        if self._args.penalty is not None:
            for line in self._args.penalty:
                line, _, comment = line.partition('#')
                line = line.strip()
                if line == '':
                    continue
                rxnid, weight = line.split(None, 1)
                weights[rxnid] = float(weight)

        # Run Fastcore and print the induced reaction set
        logger.info('Calculating Fastcore induced set on model')
        core = set(self._mm.reactions)

        induced = fastcore.fastcore(model_complete, core, epsilon,
                                    weights=weights, solver=solver)
        logger.info('Result: |A| = {}, A = {}'.format(len(induced), induced))
        added_reactions = induced - core
        logger.info('Extended: |E| = {}, E = {}'.format(
            len(added_reactions), added_reactions))

        if self._args.reaction is not None:
            maximized_reaction = self._args.reaction
        else:
            maximized_reaction = self._model.get_biomass_reaction()
            if maximized_reaction is None:
                raise ValueError('The maximized reaction was not specified')

        if not not self._mm.has_reaction(maximized_reaction):
            raise ValueError(('The biomass reaction is not a valid model' +
                              ' reaction: {}').format(maximized_reaction))

        logger.info('Flux balance on induced model maximizing {}'.format(
            maximized_reaction))
        model_induced = self._mm.copy()
        for rxnid in induced:
            model_induced.add_reaction(rxnid)
        for rxnid, flux in sorted(fluxanalysis.flux_balance(
                model_induced, maximized_reaction, tfba=enable_tfba,
                solver=solver)):
            reaction_class = 'Dbase'
            weight = weights.get(rxnid, 1)
            if self._mm.has_reaction(rxnid):
                reaction_class = 'Model'
                weight = 0
            reaction = model_complete.get_reaction(rxnid).translated_compounds(lambda x: compound_name.get(x, x))
            print('{}\t{}\t{}\t{}\t{}'.format(
                rxnid, reaction_class, weight, flux, reaction))

        logger.info('Calculating Fastcc consistent subset of induced model')
        consistent_core = fastcore.fastcc_consistent_subset(
            model_induced, epsilon, solver=solver)
        logger.info('Result: |A| = {}, A = {}'.format(
            len(consistent_core), consistent_core))
        removed_reactions = set(model_induced.reactions) - consistent_core
        logger.info('Removed: |R| = {}, R = {}'.format(
            len(removed_reactions), removed_reactions))


class FluxBalanceCommand(SolverCommandMixin, Command):
    """Run flux balance analysis on a metabolic model"""

    name = 'fba'
    title = 'Run flux balance analysis on a metabolic model'

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--no-tfba', help='Disable thermodynamic constraints on FBA',
            action='store_true')
        parser.add_argument(
            '--epsilon', type=float, help='Threshold for flux minimization',
            default=1e-5)
        parser.add_argument('reaction', help='Reaction to maximize', nargs='?')
        super(FluxBalanceCommand, cls).init_parser(parser)

    def run(self):
        """Run flux analysis command"""

        # Load compound information
        compound_name = {}
        for compound in self._model.parse_compounds():
            if 'name' in compound.properties:
                compound_name[compound.id] = compound.properties['name']
            elif compound.id not in compound_name:
                compound_name[compound.id] = compound.id

        if self._args.reaction is not None:
            reaction = self._args.reaction
        else:
            reaction = self._model.get_biomass_reaction()
            if reaction is None:
                raise ValueError('The biomass reaction was not specified')

        if not self._mm.has_reaction(reaction):
            raise ValueError('Specified reaction is not in model: {}'.format(
                reaction))

        if self._args.no_tfba:
            result = self.run_fba_minimized(reaction)
        else:
            result = self.run_tfba(reaction)

        optimum = None
        for reaction_id, fba_flux, flux in sorted(result):
            rx = self._mm.get_reaction(reaction_id)
            print('{}\t{}\t{}\t{}'.format(
                reaction_id, fba_flux, flux,
                rx.translated_compounds(lambda x: compound_name.get(x, x))))
            # Remember flux of requested reaction
            if reaction_id == reaction:
                optimum = flux

        logger.info('Maximum flux: {}'.format(optimum))

    def run_fba_minimized(self, reaction):
        """Run normal FBA and flux minimization on model, then print output"""

        solver = self._get_solver()
        fba_fluxes = dict(fluxanalysis.flux_balance(
            self._mm, reaction, tfba=False, solver=solver))
        optimum = fba_fluxes[reaction]
        epsilon = self._args.epsilon

        # Run flux minimization
        fmin_fluxes = dict(fluxanalysis.flux_minimization(
            self._mm, {reaction: optimum}, solver=solver))
        count = 0
        for reaction_id, flux in fmin_fluxes.iteritems():
            if fba_fluxes[reaction_id] - epsilon > flux:
                count += 1
            yield reaction_id, fba_fluxes[reaction_id], flux
        logger.info('Minimized reactions: {}'.format(count))

    def run_tfba(self, reaction):
        """Run FBA and tFBA on model"""

        solver = self._get_solver(integer=True)

        fba_fluxes = dict(fluxanalysis.flux_balance(
            self._mm, reaction, tfba=False, solver=solver))
        fluxes = dict(fluxanalysis.flux_balance(
            self._mm, reaction, tfba=True, solver=solver))

        for reaction_id, flux in fluxes.iteritems():
            yield reaction_id, fba_fluxes[reaction_id], flux


class FluxConsistencyCommand(SolverCommandMixin, Command):
    """Check that reactions are flux consistent in a model

    A reaction is flux consistent if there exists any steady-state flux
    solution where the flux of the given reaction is non-zero. The
    bounds on the exchange reactions can be removed when performing the
    consistency check by providing the ``--unrestricted`` option.
    """

    name = 'fluxcheck'
    title = 'Check that the model is flux consistent'

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--fastcore', help='Enable use of Fastcore algorithm',
            action='store_true')
        parser.add_argument(
            '--no-tfba',
            help='Disable thermodynamic constraints on flux check',
            action='store_true')
        parser.add_argument(
            '--epsilon', type=float, help='Flux threshold',
            default=1e-5)
        parser.add_argument(
            '--unrestricted', action='store_true',
            help='Remove limits on exchange reactions before checking'
        )
        super(FluxConsistencyCommand, cls).init_parser(parser)

    def run(self):
        """Run flux consistency check command"""

        # Load compound information
        compound_name = {}
        for compound in self._model.parse_compounds():
            if 'name' in compound.properties:
                compound_name[compound.id] = compound.properties['name']
            elif compound.id not in compound_name:
                compound_name[compound.id] = compound.id

        epsilon = self._args.epsilon

        if self._args.unrestricted:
            # Allow all exchange reactions with no flux limits
            for reaction in self._mm.reactions:
                if self._mm.is_exchange(reaction):
                    del self._mm.limits[reaction].bounds

        if self._args.fastcore:
            solver = self._get_solver()
            inconsistent = set(fastcore.fastcc(
                self._mm, epsilon, solver=solver))
        else:
            enable_tfba = not self._args.no_tfba
            if enable_tfba:
                solver = self._get_solver(integer=True)
            else:
                solver = self._get_solver()
            inconsistent = set(
                fluxanalysis.consistency_check(
                    self._mm, self._mm.reactions, epsilon,
                    tfba=enable_tfba, solver=solver))

        # Count the number of reactions that are fixed at zero. While these
        # reactions are still inconsistent, they are inconsistent because they
        # have been explicitly disabled.
        disabled_exchange = 0
        disabled_internal = 0

        count_exchange = 0
        total_exchange = 0

        count_internal = 0
        total_internal = 0

        # Print result
        for reaction in sorted(self._mm.reactions):
            disabled = self._mm.limits[reaction].bounds == (0, 0)

            if self._mm.is_exchange(reaction):
                total_exchange += 1
                count_exchange += int(reaction in inconsistent)
                disabled_exchange += int(disabled)
            else:
                total_internal += 1
                count_internal += int(reaction in inconsistent)
                disabled_internal += int(disabled)

            if reaction in inconsistent:
                rx = self._mm.get_reaction(reaction)
                rxt = rx.translated_compounds(lambda x: compound_name.get(x, x))
                print('{}\t{}'.format(reaction, rxt))

        logger.info('Model has {}/{} inconsistent internal reactions'
                    ' ({} disabled by user)'.format(
            count_internal, total_internal, disabled_internal))
        logger.info('Model has {}/{} inconsistent exchange reactions'
                    ' ({} disabled by user)'.format(
            count_exchange, total_exchange, disabled_exchange))


class FluxVariabilityCommand(SolverCommandMixin, Command):
    """Run flux variablity analysis on a metabolic model"""

    name = 'fva'
    title = 'Run flux variability analysis on a metabolic model'

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--no-tfba', help='Disable thermodynamic constraints on FVA',
            action='store_true')
        parser.add_argument('reaction', help='Reaction to maximize', nargs='?')
        super(FluxVariabilityCommand, cls).init_parser(parser)

    def run(self):
        """Run flux variability command"""

        # Load compound information
        compound_name = {}
        for compound in self._model.parse_compounds():
            if 'name' in compound.properties:
                compound_name[compound.id] = compound.properties['name']
            elif compound.id not in compound_name:
                compound_name[compound.id] = compound.id

        if self._args.reaction is not None:
            reaction = self._args.reaction
        else:
            reaction = self._model.get_biomass_reaction()
            if reaction is None:
                raise ValueError('The biomass reaction was not specified')

        if not self._mm.has_reaction(reaction):
            raise ValueError('Specified reaction is not in model: {}'.format(
                reaction))

        enable_tfba = not self._args.no_tfba
        if enable_tfba:
            solver = self._get_solver(integer=True)
        else:
            solver = self._get_solver()

        fba_fluxes = dict(fluxanalysis.flux_balance(
            self._mm, reaction, tfba=False, solver=solver))
        optimum = fba_fluxes[reaction]

        flux_bounds = fluxanalysis.flux_variability(
            self._mm, sorted(self._mm.reactions), {reaction: optimum},
            tfba=enable_tfba, solver=solver)
        for reaction_id, bounds in flux_bounds:
            rx = self._mm.get_reaction(reaction_id)
            rxt = rx.translated_compounds(lambda x: compound_name.get(x, x))
            print('{}\t{}\t{}\t{}'.format(
                reaction_id, bounds[0], bounds[1], rxt))


class FormulaBalanceCommand(Command):
    """Check whether reactions in a given database or model are balanced

    Balanced reactions are those reactions where the number of atoms
    is consistent on the left and right side of the reaction equation.
    Reactions that are not balanced will be printed out.
    """

    name = 'formulacheck'
    title = 'Check formula balance on a model or database'

    def run(self):
        """Run formula balance command"""

        # Mapping from compound id to formula
        compound_formula = {}
        for compound in self._model.parse_compounds():
            # Only cpd11632 (Photon) is allowed to have an empty formula.
            if compound.formula is None:
                if compound.id == 'cpd11632':
                    compound_formula[compound.id] = Formula()
            else:
                try:
                    f = Formula.parse(compound.formula).flattened()
                    compound_formula[compound.id] = f
                except ValueError as e:
                    logger.warning(
                        'Error parsing {}: {}'.format(compound.formula, e))

        # Create a set of known mass-inconsistent reactions
        exchange = set()
        for reaction_id in self._mm.reactions:
            if self._mm.is_exchange(reaction_id):
                exchange.add(reaction_id)

        def multiply_formula(compound_list):
            for compound, count in compound_list:
                yield count * compound_formula.get(compound.name, Formula())

        count = 0
        unbalanced = 0
        unchecked = 0
        for reaction in self._mm.reactions:
            if reaction in exchange:
                continue

            count += 1
            rx = self._mm.get_reaction(reaction)

            # Skip reaction if any compounds have undefined formula
            for compound, _ in rx.compounds:
                if compound.name not in compound_formula:
                    unchecked += 1
                    break
            else:
                left_form = reduce(
                    operator.or_, multiply_formula(rx.left), Formula())
                right_form = reduce(
                    operator.or_, multiply_formula(rx.right), Formula())

                if right_form != left_form:
                    unbalanced += 1
                    right_missing, left_missing = Formula.balance(
                        right_form, left_form)

                    print('{}\t{}\t{}\t{}\t{}'.format(
                        reaction, left_form, right_form,
                        left_missing, right_missing))

        logger.info('Unbalanced reactions: {}/{}'.format(unbalanced, count))
        logger.info('Unchecked reactions due to missing formula: {}/{}'.format(
            unchecked, count))


class GapFillCommand(SolverCommandMixin, Command):
    """Command that runs GapFind and GapFill on a metabolic model"""

    name = 'gapfill'
    title = 'Run GapFind and GapFill on a metabolic model'

    def run(self):
        """Run GapFill command"""

        # Load compound information
        compound_name = {}
        for compound in self._model.parse_compounds():
            if 'name' in compound.properties:
                compound_name[compound.id] = compound.properties['name']
            elif compound.id not in compound_name:
                compound_name[compound.id] = compound.id

        model_compartments = set(self._mm.compartments)

        solver = self._get_solver(integer=True)

        # Run GapFind on model
        logger.info('Searching for blocked compounds')
        blocked = set(compound for compound in gapfind(self._mm, solver=solver)
                      if compound.compartment is not 'e')
        if len(blocked) > 0:
            logger.info('Blocked compounds')
            for compound in blocked:
                print(compound.translate(lambda x: compound_name.get(x, x)))

        if len(blocked) > 0:
            # Add exchange and transport reactions to database
            model_complete = self._mm.copy()
            logger.info('Adding database, exchange and transport reactions')
            model_complete.add_all_database_reactions(model_compartments)
            model_complete.add_all_exchange_reactions()
            model_complete.add_all_transport_reactions()

            logger.info('Searching for reactions to fill gaps')
            added_reactions, reversed_reactions = gapfill(
                model_complete, self._mm.reactions, blocked, solver=solver)

            for rxnid in added_reactions:
                rx = model_complete.get_reaction(rxnid)
                rxt = rx.translated_compounds(
                    lambda x: compound_name.get(x, x))
                print('{}\t{}\t{}'.format(rxnid, 'Add', rxt))

            for rxnid in reversed_reactions:
                rx = model_complete.get_reaction(rxnid)
                rxt = rx.translated_compounds(
                    lambda x: compound_name.get(x, x))
                print('{}\t{}\t{}'.format(rxnid, 'Reverse', rxt))
        else:
            logger.info('No blocked compounds found')


class MassConsistencyCommand(SolverCommandMixin, Command):
    """Command that checks whether a database is mass consistent"""

    name = 'masscheck'
    title = 'Run mass consistency check on a database'

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--exclude', metavar='reaction', action='append', type=str,
            default=[], help='Exclude reaction from mass consistency')
        parser.add_argument(
            '--epsilon', type=float, help='Mass threshold',
            default=1e-5)
        parser.add_argument(
            '--checked', metavar='reaction', action='append', type=str,
            default=[], help='Mark reaction as already checked (no residual)')
        super(MassConsistencyCommand, cls).init_parser(parser)

    def run(self):
        # Load compound information
        compound_name = {}
        zeromass = set()
        for compound in self._model.parse_compounds():
            if 'name' in compound.properties:
                compound_name[compound.id] = compound.properties['name']
            elif compound.id not in compound_name:
                compound_name[compound.id] = compound.id

            if compound.properties.get('zeromass', False):
                zeromass.add(compound.id)

        # Create a set of known mass-inconsistent reactions
        exchange = set()
        for reaction_id in self._mm.reactions:
            if self._mm.is_exchange(reaction_id):
                exchange.add(reaction_id)

        # Other reactions to exclude from consistency check
        exclude = set(self._args.exclude)

        # Add biomass reaction to be excluded
        biomass_reaction = self._model.get_biomass_reaction()
        if biomass_reaction is not None:
            exclude.add(biomass_reaction)

        # Create set of checked reactions
        checked = set(self._args.checked)

        solver = self._get_solver()

        known_inconsistent = exclude | exchange

        logger.info('Mass consistency on database')
        epsilon = self._args.epsilon
        compound_iter = massconsistency.check_compound_consistency(
            self._mm, solver, known_inconsistent, zeromass)

        good = 0
        total = 0
        for compound, mass in sorted(compound_iter, key=lambda x: (x[1], x[0]),
                                     reverse=True):
            if mass >= 1-epsilon or compound.name in zeromass:
                good += 1
            total += 1
            print('{}: {}'.format(compound.translate(
                lambda x: compound_name.get(x, x)), mass))
        logger.info('Consistent compounds: {}/{}'.format(good, total))

        logger.info('Is consistent? {}'.format(
            massconsistency.is_consistent(
                self._mm, solver, known_inconsistent, zeromass)))

        good = 0
        total = 0
        reaction_iter, compound_iter = (
            massconsistency.check_reaction_consistency(
                self._mm, solver=solver, exchange=known_inconsistent,
                checked=checked, zeromass=zeromass))
        for reaction_id, residual in sorted(
                reaction_iter, key=lambda x: abs(x[1]), reverse=True):
            total += 1
            if abs(residual) >= epsilon:
                reaction = self._mm.get_reaction(reaction_id)
                rxt = reaction.translated_compounds(
                    lambda x: compound_name.get(x, x))
                print('{}\t{}\t{}'.format(reaction_id, residual, rxt))
            else:
                good += 1
        logger.info('Consistent reactions: {}/{}'.format(good, total))


class RandomSparseNetworkCommand(SolverCommandMixin, Command):
    """Find random minimal network of model

    Given a reaction to optimize and a threshold, delete reactions randomly
    until the flux of the reaction to optimize falls under the threshold.
    Keep deleting reactions until no more reactions can be deleted. By default
    this uses standard FBA (not tFBA). Since the internal fluxes are irrelevant
    the FBA and tFBA are equivalent for this purpose.
    """

    name = 'randomsparse'
    title = 'Generate a random sparse network of model reactions'

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--tfba', help='Enable thermodynamic constraints on FBA',
            action='store_true', default=False)
        parser.add_argument(
            '--reaction', help='Reaction to maximize', nargs='?')
        parser.add_argument(
            'threshold', help='Relative threshold of max reaction flux',
            type=float)
        parser.add_argument(
            '--exchange', help='Only analyze the exchange reaction in model',
            action='store_true')
        super(RandomSparseNetworkCommand, cls).init_parser(parser)

    def run(self):
        if self._args.reaction is not None:
            reaction = self._args.reaction
        else:
            reaction = self._model.get_biomass_reaction()
            if reaction is None:
                raise ValueError('The biomass reaction was not specified')

        if not self._mm.has_reaction(reaction):
            raise ValueError('Specified reaction is not in model: {}'.format(
                reaction))

        threshold = self._args.threshold
        if threshold < 0.0 or threshold > 1.0:
            raise ValueError(
                'Invalid threshold, must be in [0;1]: {}'.format(threshold))

        if not self._args.tfba:
            fb_problem = fluxanalysis.FluxBalanceProblem
            solver = self._get_solver()
        else:
            fb_problem = fluxanalysis.FluxBalanceTDProblem
            solver = self._get_solver(integer=True)

        p = fb_problem(self._mm, solver)
        p.solve(reaction)
        flux_threshold = p.get_flux(reaction) * threshold

        logger.info('Flux threshold for {} is {}'.format(reaction, flux_threshold))

        if self._args.exchange:
            model_test = self._mm.copy()
            essential = { reaction }
            deleted = set()
            exchange = set()
            for reaction_id in self._mm.reactions:
                if self._mm.is_exchange(reaction_id):
                    exchange.add(reaction_id)
            test_set = set(exchange) - essential
        else:
            model_test = self._mm.copy()
            essential = { reaction }
            deleted = set()
            test_set = set(self._mm.reactions) - essential

        while len(test_set) > 0:
            testing_reaction = random.sample(test_set, 1)[0]
            test_set.remove(testing_reaction)
            saved_bounds = model_test.limits[testing_reaction].bounds
            model_test.limits[testing_reaction].bounds = 0, 0

            logger.info('Trying FBA without reaction {}...'.format(testing_reaction))

            p = fb_problem(model_test, solver)
            try:
                p.solve(reaction)
            except fluxanalysis.FluxBalanceError:
                logger.info('FBA is infeasible, marking {} as essential'.format(testing_reaction))
                model_test.limits[testing_reaction].bound = saved_bounds
                essential.add(testing_reaction)
                continue

            logger.debug('Reaction {} has flux {}'.format(reaction, p.get_flux(reaction)))

            if p.get_flux(reaction) < flux_threshold:
                model_test.limits[testing_reaction].bounds = saved_bounds
                essential.add(testing_reaction)
                logger.info('Reaction {} was essential'.format(testing_reaction))
            else:
                deleted.add(testing_reaction)
                logger.info('Reaction {} was deleted'.format(testing_reaction))
        if self._args.exchange:
            for reaction_id in sorted(exchange):
                value = 0 if reaction_id in deleted else 1
                print('{}\t{}'.format(reaction_id, value))
        else:
            for reaction_id in sorted(self._mm.reactions):
                value = 0 if reaction_id in deleted else 1
                print('{}\t{}'.format(reaction_id, value))

class RobustnessCommand(SolverCommandMixin, Command):
    """Run robustness analysis on metabolic model

    Given a reaction to maximize and a reaction to vary,
    the robustness analysis will run FBA while fixing the
    reaction to vary at each iteration. The reaction will
    be fixed at the specified number of steps between the
    minimum and maximum flux value specified in the model.
    """

    name = 'robustness'
    title = 'Run robustness analysis on a metabolic model'

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--steps', metavar='N', type=int, default=10,
            help='Number of flux value steps for varying reaction')
        parser.add_argument(
            '--minimum', metavar='V', type=float,
            help='Minumum flux value of varying reacton')
        parser.add_argument(
            '--maximum', metavar='V', type=float,
            help='Maximum flux value of varying reacton')
        parser.add_argument(
            '--no-tfba', help='Disable thermodynamic constraints on FBA',
            action='store_true')
        parser.add_argument(
            '--reaction', help='Reaction to maximize', nargs='?')
        parser.add_argument('varying', help='Reaction to vary')
        super(RobustnessCommand, cls).init_parser(parser)

    def run(self):
        """Run flux analysis command"""

        def run_fba_fmin(model, reaction):
            fba = fluxanalysis.FluxBalanceProblem(model, solver=solver)
            fba.solve(reaction)
            optimum = fba.get_flux(reaction)
            return fluxanalysis.flux_minimization(
                model, {reaction: optimum}, solver=solver)

        def run_tfba(model, reaction):
            return fluxanalysis.flux_balance(
                model, reaction, tfba=True, solver=solver)

        if self._args.no_tfba:
            run_fba = run_fba_fmin
            solver = self._get_solver()
        else:
            run_fba = run_tfba
            solver = self._get_solver(integer=True)

        # Load compound information
        compound_name = {}
        for compound in self._model.parse_compounds():
            if 'name' in compound.properties:
                compound_name[compound.id] = compound.properties['name']
            elif compound.id not in compound_name:
                compound_name[compound.id] = compound.id

        if self._args.reaction is not None:
            reaction = self._args.reaction
        else:
            reaction = self._model.get_biomass_reaction()
            if reaction is None:
                raise ValueError('The biomass reaction was not specified')

        if not self._mm.has_reaction(reaction):
            raise ValueError('Specified reaction is not in model: {}'.format(
                reaction))

        varying_reaction = self._args.varying
        if not self._mm.has_reaction(varying_reaction):
            raise ValueError('Specified reaction is not in model: {}'.format(
                varying_reaction))

        steps = self._args.steps
        if steps <= 0:
            raise ValueError('Invalid number of steps: {}\n'.format(steps))

        # Run FBA on model at different fixed flux values
        flux_min = self._mm.limits[varying_reaction].lower
        flux_max = self._mm.limits[varying_reaction].upper
        if self._args.minimum is not None:
            flux_min = self._args.minimum
        if self._args.maximum is not None:
            flux_max = self._args.maximum

        if flux_min > flux_max:
            raise ValueError('Invalid flux range: {}, {}\n'.format(
                flux_min, flux_max))

        for i in xrange(steps):
            fixed_flux = flux_min + i*(flux_max - flux_min)/float(steps-1)
            test_model = self._mm.copy()
            test_model.limits[varying_reaction].bounds = fixed_flux, fixed_flux

            try:
                for other_reaction, flux in run_fba(test_model, reaction):
                    print('{}\t{}\t{}'.format(
                        other_reaction, fixed_flux, flux))
            except fluxanalysis.FluxBalanceError:
                pass


class SBMLExport(Command):
    """Export model as SBML file"""

    name = 'sbmlexport'
    title = 'Export model as SBML file'

    def run(self):
        writer = sbml.SBMLWriter()
        writer.write_model(
            sys.stdout, self._mm, self._model.parse_compounds())


class SearchCommand(Command):
    """Defines the search command"""

    name = 'search'
    title = 'Search the database of reactions or compounds'

    @classmethod
    def init_parser(cls, parser):
        """Initialize argument parser"""
        subparsers = parser.add_subparsers(title='Search domain')

        # Compound subcommand
        parser_compound = subparsers.add_parser(
            'compound', help='Search in compounds')
        parser_compound.set_defaults(which='compound')
        parser_compound.add_argument(
            '--id', '-i', dest='id', metavar='id', type=str,
            default=[], action='append', help='Compound ID')
        parser_compound.add_argument(
            '--name', '-n', dest='name', metavar='name', type=str,
            default=[], action='append', help='Name of compound')

        # Reaction subcommand
        parser_reaction = subparsers.add_parser(
            'reaction', help='Search in reactions')
        parser_reaction.set_defaults(which='reaction')
        parser_reaction.add_argument(
            '--id', '-i', dest='id', metavar='id', type=str,
            default=[], action='append', help='Reaction ID')
        parser_reaction.add_argument(
            '--compound', '-c', dest='compound', metavar='compound', type=str,
            default=[], action='append',
            help='Comma-separated list of compound IDs')

    def run(self):
        """Run search command"""

        def parse_compound(s):
            """Parse a compound specification with optional compartment"""
            m = re.match(r'([^\[]+)\[(\w+)\]', s)
            if m is not None:
                return Compound(m.group(1), compartment=m.group(2))
            return Compound(s)

        def filter_search_term(s):
            return re.sub(r'[^a-z0-9]+', '', s.lower())

        # Load compound information
        compound_name = {}
        compound_synonyms = {}
        compound_formula = {}
        compound_for_name = {}
        for compound in self._model.parse_compounds():
            if 'name' in compound.properties:
                compound_name[compound.id] = compound.properties['name']
            elif compound.id not in compound_name:
                compound_name[compound.id] = compound.id

            compound_synonyms.setdefault(compound.id, []).extend(
                compound.properties.get('names', []))

            if compound.formula is not None and not '.' in compound.formula:
                compound_formula[compound.id] = Formula.parse(compound.formula)

        # Create references from names to id
        for compound_id in compound_name.iterkeys():
            names = ([compound_name[compound_id]] +
                compound_synonyms[compound_id])

            for n in names:
                n = filter_search_term(n)
                compound_for_name[n] = compound_id

        which_command = self._args.which
        if which_command == 'compound':
            compound_ids = []
            for name in self._args.name:
                search_name = filter_search_term(name)
                if search_name in compound_for_name:
                    compound_ids.append(compound_for_name[search_name])
            for compound_id in self._args.id:
                if compound_id in compound_name:
                    compound_ids.append(compound_id)

            # Show results
            for compound_id in compound_ids:
                print('ID: {}'.format(compound_id))
                if compound_id in compound_name:
                    print('Name: {}'.format(compound_name[compound_id]))
                if compound_id in compound_synonyms:
                    print('Synonyms: {}'.format(
                        ', '.join(compound_synonyms[compound_id])))
                if compound_id in compound_formula:
                    print('Formula: {}'.format(compound_formula[compound_id]))
                print()
        elif which_command == 'reaction':
            reaction_ids = []
            for reaction_id in self._args.id:
                if reaction_id in self._mm.reactions:
                    reaction_ids.append(reaction_id)
            for compound_list in self._args.compound:
                matched_reactions = None
                for compound_spec in compound_list.split(','):
                    compound = parse_compound(compound_spec.strip())
                    if matched_reactions is None:
                        matched_reactions = set(
                            self._mm.get_compound_reactions(compound))
                    else:
                        matched_reactions &= set(
                            self._mm.get_compound_reactions(compound))
                if matched_reactions is not None:
                    reaction_ids.extend(matched_reactions)

            # Show results
            for reaction_id in reaction_ids:
                print('ID: {}'.format(reaction_id))

                reaction = self._mm.get_reaction(reaction_id)
                print('Reaction (IDs): {}'.format(
                    self._mm.get_reaction(reaction_id)))
                translated_reaction = reaction.translated_compounds(
                    lambda x: compound_name.get(x, x))
                if reaction != translated_reaction:
                    print('Reaction (names): {}'.format(translated_reaction))
                print()


def main(command_class=None):
    """Run the command line interface with the given :class:`Command`

    If no command class is specified the user will be able to select a specific
    command through the first command line argument.
    """

    # Set up logging for the command line interface
    if 'DEBUG' in os.environ:
        level = getattr(logging, os.environ['DEBUG'].upper(), None)
        if level is not None:
            logging.basicConfig(level=level)
    else:
        logging.basicConfig(level=logging.INFO)

    title = 'Metabolic modeling tools'
    if command_class is not None:
        title = command_class.title

    parser = argparse.ArgumentParser(description=title)
    parser.add_argument('--model', metavar='file', default='.',
                        help='Model definition')
    parser.add_argument(
        '-V', '--version', action='version',
        version='%(prog)s ' + package_version)

    if command_class is not None:
        # Command explicitly given, only allow that command
        command_class.init_parser(parser)
        parser.set_defaults(command=command_class)
    else:
        # Discover all available commands
        commands = {}
        for command_class in Command.__subclasses__():
            command_name = getattr(command_class, 'name', None)
            command_title = getattr(command_class, 'title', None)
            if (command_name is not None and
                    command_title is not None and
                    command_name not in commands):
                commands[command_name] = (command_class.title,
                                          command_class)

        # Create parsers for subcommands
        subparsers = parser.add_subparsers(title='Command')
        for name, values in sorted(commands.iteritems()):
            title, command_class = values
            subparser = subparsers.add_parser(name, help=title)
            subparser.set_defaults(command=command_class)
            command_class.init_parser(subparser)

    args = parser.parse_args()

    # Load model definition
    model = NativeModel(args.model)

    # Instantiate command with model and run
    command = args.command(model, args)
    command.run()


if __name__ == '__main__':
    main()
