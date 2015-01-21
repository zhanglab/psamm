
'''Utilities for the command line interface'''

import sys
import argparse
import operator
import re
import logging

from .fastcore import Fastcore
from .formula import Formula, Radical
from .gapfill import gapfind, gapfill
from .massconsistency import MassConsistencyCheck
from .database import DictDatabase, ChainedDatabase
from .metabolicmodel import MetabolicModel
from .reaction import Reaction, Compound
from .datasource import internal, modelseed
from . import fluxanalysis

# Module-level logging
logger = logging.getLogger(__name__)

class Command(object):
    '''Represents a command in the interface, operating on a model or database

    Subclasses must define name and title as class attributes, and implement
    __call__ to handle command execution. Arguments will be given as keyword
    arguments. The keywords database, model and compounds will be set.

    In addition, init_parser() can be implemented which will allow the
    command to initialize an instance of ArgumentParser as desired. The
    resulting arguments will be given as keyword arguments to __call__.
    '''

    def init_parser(self, parser):
        '''Initialize command line parser (argparse.ArgumentParser)'''
        pass

    def __call__(self):
        '''Execute command'''
        pass

class ChargeBalanceCommand(Command):
    '''Check whether compound charge in a given database or model is balanced

    Balanced reactions are those reactions where the total charge
    is consistent on the left and right side of the reaction equation.
    Reactions that are not balanced will be printed out.'''

    name = 'chargecheck'
    title = 'Check charge balance on a model or database'

    def init_parser(self, parser):
        parser.add_argument('charge', metavar='file', type=argparse.FileType('r'),
                            help='List of charges for database compounds')

    def __call__(self, model, compounds, **kwargs):
        '''Run charge balance command'''

        # Load compound information
        compound_name = {}
        for compound in compounds:
            compound_name[compound.id] = compound.name if compound.name is not None else compound.id

        # Mapping from compound id to charge
        compound_charge = {}
        for line in kwargs['charge']:
            line, _, comment = line.partition('#')
            line = line.strip()
            if line == '':
                continue
            compound_id, charge = line.split(None, 1)
            compound_charge[compound_id] = int(charge)

        # Create a set of known charge-inconsistent reactions
        exchange = set()
        for reaction_id in model.reactions:
            if model.is_exchange(reaction_id):
                exchange.add(reaction_id)

        def reaction_charges(reaction_id):
            for compound, value in model.get_reaction_values(reaction_id):
                yield compound_charge.get(compound.name, 0) * value

        count = 0
        for reaction in sorted(model.reactions):
            if reaction not in exchange:
                charge = sum(reaction_charges(reaction))
                if charge != 0:
                    count += 1
                    rx = model.get_reaction(reaction).translated_compounds(lambda x: compound_name.get(x, x))
                    print '{}\t{}\t{}'.format(reaction, charge, rx)

        logger.info('Unbalanced reactions: {}'.format(count))

class ConsoleCommand(Command):
    '''Start an interactive Python console with the given model loaded'''

    name = 'console'
    title = 'Start Python console with metabolic model loaded'

    def init_parser(self, parser):
        parser.add_argument('--type', choices=('python', 'ipython', 'ipython-kernel'),
                            default='python', help='type of console to open')

    def open_python(self, message, namespace):
        '''Open interactive python console'''

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

    def __call__(self, **kwargs):
        message = 'Metabolic model has been loaded into "model" and "compounds".'
        console_type = kwargs['type']

        if console_type == 'python':
            self.open_python(message, kwargs)
        elif console_type == 'ipython':
            self.open_ipython(message, kwargs)
        elif console_type == 'ipython-kernel':
            self.open_ipython_kernel(message, kwargs)

class FastGapFillCommand(Command):
    '''Run FastGapFill algorithm on a metabolic model'''

    name = 'fastgapfill'
    title = 'Run FastGapFill on a metabolic model'

    def init_parser(self, parser):
        parser.add_argument('--penalty', metavar='file', type=argparse.FileType('r'),
                            help='List of penalty scores for database reactions')
        parser.add_argument('reaction', help='Reaction to maximize')

    def __call__(self, model, compounds, **kwargs):
        '''Run FastGapFill command'''

        # Create fastcore object
        from metnet.lpsolver import cplex
        solver = cplex.Solver()
        fastcore = Fastcore(solver)

        # Load compound information
        compound_name = {}
        for compound in compounds:
            compound_name[compound.id] = compound.name if compound.name is not None else compound.id

        epsilon = 1e-5
        model_compartments = { None, 'e' }

        # Add exchange and transport reactions to database
        model_complete = model.copy()
        logger.info('Adding database, exchange and transport reactions')
        db_added = model_complete.add_all_database_reactions(model_compartments)
        ex_added = model_complete.add_all_exchange_reactions()
        tp_added = model_complete.add_all_transport_reactions()

        # Add penalty weights on reactions
        weights = {}
        #weights.update((rxnid, 25) for rxnid in db_added)
        weights.update((rxnid, 50) for rxnid in tp_added)
        weights.update((rxnid, 250) for rxnid in ex_added)

        if kwargs.get('penalty') is not None:
            for line in kwargs['penalty']:
                line, _, comment = line.partition('#')
                line = line.strip()
                if line == '':
                    continue
                rxnid, weight = line.split(None, 1)
                weights[rxnid] = float(weight)

        # Run Fastcore and print the induced reaction set
        logger.info('Calculating Fastcore induced set on model')
        core = set(model.reactions)

        induced = fastcore.fastcore(model_complete, core, epsilon, weights=weights)
        logger.info('Result: |A| = {}, A = {}'.format(len(induced), induced))
        added_reactions = induced - core
        logger.info('Extended: |E| = {}, E = {}'.format(len(added_reactions), added_reactions))

        maximized_reaction = kwargs['reaction']
        logger.info('Flux balance on induced model maximizing {}'.format(maximized_reaction))
        model_induced = model.copy()
        for rxnid in induced:
            model_induced.add_reaction(rxnid)
        for rxnid, flux in sorted(fluxanalysis.flux_balance(model_induced, maximized_reaction, solver=solver)):
            reaction_class = 'Dbase'
            weight = weights.get(rxnid, 1)
            if model.has_reaction(rxnid):
                reaction_class = 'Model'
                weight = 0
            reaction = model_complete.get_reaction(rxnid).translated_compounds(lambda x: compound_name.get(x, x))
            print '{}\t{}\t{}\t{}\t{}'.format(rxnid, reaction_class, weight, flux, reaction)

        logger.info('Calculating Fastcc consistent subset of induced model')
        consistent_core = fastcore.fastcc_consistent_subset(model_induced, epsilon)
        logger.info('Result: |A| = {}, A = {}'.format(len(consistent_core), consistent_core))
        removed_reactions = set(model_induced.reactions) - consistent_core
        logger.info('Removed: |R| = {}, R = {}'.format(len(removed_reactions), removed_reactions))

class FluxBalanceCommand(Command):
    '''Run flux balance analysis on a metabolic model'''

    name = 'fba'
    title = 'Run flux balance analysis on a metabolic model'

    def init_parser(self, parser):
        parser.add_argument('--no-tfba', help='Disable thermodynamic constraints on FBA',
                            action='store_true')
        parser.add_argument('reaction', help='Reaction to maximize')

    def __call__(self, model, compounds, **kwargs):
        '''Run flux analysis command'''

        # Load compound information
        compound_name = {}
        for compound in compounds:
            compound_name[compound.id] = compound.name if compound.name is not None else compound.id

        reaction = kwargs['reaction']
        if not model.has_reaction(reaction):
            raise ValueError('Specified reaction is not in model: {}'.format(reaction))

        if kwargs['no_tfba']:
            result = self.run_fba_minimized(model, reaction)
        else:
            result = self.run_tfba(model, reaction)

        optimum = None
        for reaction_id, fba_flux, flux in sorted(result):
            rx = model.get_reaction(reaction_id)
            print '{}\t{}\t{}\t{}'.format(reaction_id, fba_flux, flux,
                                            rx.translated_compounds(lambda x: compound_name.get(x, x)))
            # Remember flux of requested reaction
            if reaction_id == reaction:
                optimum = flux

        logger.info('Maximum flux: {}'.format(optimum))

    def run_fba_minimized(self, model, reaction):
        '''Run normal FBA and flux minimization on model, then print output'''

        from .lpsolver import cplex
        solver = cplex.Solver()

        fba_fluxes = dict(fluxanalysis.flux_balance(model, reaction, solver=solver))
        optimum = fba_fluxes[reaction]

        epsilon = 1e-5

        # Run flux minimization
        fmin_fluxes = dict(fluxanalysis.flux_minimization(model, { reaction: optimum }, solver=solver))
        count = 0
        for reaction_id, flux in fmin_fluxes.iteritems():
            if fba_fluxes[reaction_id] - epsilon > flux:
                count += 1
            yield reaction_id, fba_fluxes[reaction_id], flux
        logger.info('Minimized reactions: {}'.format(count))

    def run_tfba(self, model, reaction):
        '''Run FBA and tFBA on model'''

        from .lpsolver import cplex
        solver = cplex.Solver()

        fba_fluxes = dict(fluxanalysis.flux_balance(model, reaction, solver=solver))
        fluxes = dict(fluxanalysis.flux_balance_td(model, reaction, solver=solver))
        for reaction_id, flux in fluxes.iteritems():
            yield reaction_id, fba_fluxes[reaction_id], flux

class FluxConsistencyCommand(Command):
    '''Check that reactions are flux consistent in a model

    A reaction is flux consistent if there exists any steady-state
    flux solution where the flux of the given reaction is non-zero.'''

    name = 'fluxconsistency'
    title = 'Check that the model is flux consistent'

    def init_parser(self, parser):
        parser.add_argument('--no-fastcore', help='Disable use of Fastcore algorithm',
                            action='store_true')

    def __call__(self, model, compounds, **kwargs):
        '''Run flux consistency check command'''

        # Load compound information
        compound_name = {}
        for compound in compounds:
            compound_name[compound.id] = compound.name if compound.name is not None else compound.id

        from metnet.lpsolver import cplex
        solver = cplex.Solver()
        epsilon = 1e-5

        if kwargs['no_fastcore']:
            inconsistent = set(fluxanalysis.consistency_check(model, model.reactions, epsilon, solver=solver))
        else:
            # Create fastcore object
            fastcore = Fastcore(solver)
            inconsistent = set(fastcore.fastcc(model, epsilon))

        # Print result
        for reaction in sorted(inconsistent):
            rx = model.get_reaction(reaction).translated_compounds(lambda x: compound_name.get(x, x))
            print '{}\t{}'.format(reaction, rx)

        logger.info('Model has {} inconsistent reactions.'.format(len(inconsistent)))

class FluxVariabilityCommand(Command):
    '''Run flux variablity analysis on a metabolic model'''

    name = 'fva'
    title = 'Run flux variability analysis on a metabolic model'

    def init_parser(self, parser):
        parser.add_argument('reaction', help='Reaction to maximize')

    def __call__(self, model, compounds, **kwargs):
        '''Run flux variability command'''

        # Load compound information
        compound_name = {}
        for compound in compounds:
            compound_name[compound.id] = compound.name if compound.name is not None else compound.id

        reaction = kwargs['reaction']
        if not model.has_reaction(reaction):
            raise ValueError('Specified reaction is not in model: {}'.format(reaction))

        from .lpsolver import cplex
        solver = cplex.Solver()

        fba_fluxes = dict(fluxanalysis.flux_balance(model, reaction, solver=solver))
        optimum = fba_fluxes[reaction]

        flux_bounds = sorted(fluxanalysis.flux_variability(model, model.reactions, { reaction: optimum},
                                                            solver=solver))
        for reaction_id, bounds in flux_bounds:
            rx = model.get_reaction(reaction_id).translated_compounds(lambda x: compound_name.get(x, x))
            print '{}\t{}\t{}\t{}'.format(reaction_id, bounds[0], bounds[1], rx)

class FormulaBalanceCommand(Command):
    '''Check whether reactions in a given database or model are balanced

    Balanced reactions are those reactions where the number of atoms
    is consistent on the left and right side of the reaction equation.
    Reactions that are not balanced will be printed out.'''

    name = 'formulacheck'
    title = 'Check formula balance on a model or database'

    def __call__(self, model, compounds):
        '''Run formula balance command'''

        # Mapping from compound id to formula
        compound_formula = {}
        for compound in compounds:
            # Create pseudo-radical group for compounds with
            # missing formula, so they don't match up. Only
            # cpd11632 (Photon) is allowed to have an empty formula.
            if compound.formula is None or '.' in compound.formula:
                if compound.id != 'cpd11632':
                    f = Formula({Radical('R'+compound.id): 1})
                else:
                    f = Formula()
            else:
                try:
                    f = Formula.parse(compound.formula).flattened()
                except ValueError as e:
                    logger.warning('Error parsing {}: {}'.format(compound.formula, e))
            compound_formula[compound.id] = f

        # Create a set of known mass-inconsistent reactions
        exchange = set()
        for reaction_id in model.reactions:
            if model.is_exchange(reaction_id):
                exchange.add(reaction_id)

        def multiply_formula(compound_list):
            for compound, count in compound_list:
                yield count * compound_formula.get(compound.name, Formula())

        for reaction in model.reactions:
            if reaction not in exchange:
                rx = model.get_reaction(reaction)
                left_form = reduce(operator.or_, multiply_formula(rx.left), Formula())
                right_form = reduce(operator.or_, multiply_formula(rx.right), Formula())

                if right_form != left_form:
                    right_missing, left_missing = Formula.balance(right_form, left_form)
                    print '{}\t{}\t{}\t{}\t{}'.format(reaction, left_form, right_form, left_missing, right_missing)

class GapFillCommand(Command):
    '''Command that runs GapFind and GapFill on a metabolic model'''

    name = 'gapfill'
    title = 'Run GapFind and GapFill on a metabolic model'

    def __call__(self, model, compounds):
        '''Run GapFill command'''

        # Load compound information
        compound_name = {}
        for compound in compounds:
            compound_name[compound.id] = compound.name if compound.name is not None else compound.id

        model_compartments = { None, 'e' }

        from .lpsolver import cplex
        solver = cplex.Solver()

        # Run GapFind on model
        logger.info('Searching for blocked compounds')
        blocked = set(compound for compound in gapfind(model, solver=solver) if compound.compartment is not 'e')
        if len(blocked) > 0:
            logger.info('Blocked compounds')
            for compound in blocked:
                print compound.translate(lambda x: compound_name.get(x, x))

        if len(blocked) > 0:
            # Add exchange and transport reactions to database
            model_complete = model.copy()
            logger.info('Adding database, exchange and transport reactions')
            model_complete.add_all_database_reactions(model_compartments)
            model_complete.add_all_exchange_reactions()
            model_complete.add_all_transport_reactions()

            logger.info('Searching for reactions to fill gaps')
            added_reactions, reversed_reactions = gapfill(model_complete, model.reactions, blocked, solver=solver)

            for rxnid in added_reactions:
                rx = model_complete.get_reaction(rxnid).translated_compounds(lambda x: compound_name.get(x, x))
                print '{}\t{}\t{}'.format(rxnid, 'Add', rx)

            for rxnid in reversed_reactions:
                rx = model_complete.get_reaction(rxnid).translated_compounds(lambda x: compound_name.get(x, x))
                print '{}\t{}\t{}'.format(rxnid, 'Reverse', rx)
        else:
            logger.info('No blocked compounds found')

class MassConsistencyCommand(Command):
    '''Command that checks whether a database is mass consistent'''

    name = 'masscheck'
    title = 'Run mass consistency check on a database'

    def init_parser(self, parser):
        parser.add_argument('--exclude', metavar='reaction', action='append',
                            type=str, default=[], help='Exclude reaction from mass consistency')

    def __call__(self, model, compounds, **kwargs):
        # Load compound information
        compound_name = {}
        for compound in compounds:
            compound_name[compound.id] = compound.name if compound.name is not None else compound.id

        # Create a set of known mass-inconsistent reactions
        exchange = set()
        for reaction_id in model.reactions:
            if model.is_exchange(reaction_id):
                exchange.add(reaction_id)

        # Other reactions to exclude from consistency check
        exclude = set(kwargs['exclude'])

        # Create set of compounds allowed to have mass zero
        zeromass = set()
        zeromass.add('cpd11632') # Photon
        zeromass.add('cpd12713') # Electron

        from .lpsolver import cplex
        solver = cplex.Solver()

        mass_consistency = MassConsistencyCheck(solver=solver)
        known_inconsistent = exclude | exchange

        logger.info('Mass consistency on database')
        epsilon = 1e-5
        compound_iter = mass_consistency.check_compound_consistency(model, known_inconsistent, zeromass)

        logger.info('Compound consistency')
        good = 0
        total = 0
        for compound, mass in sorted(compound_iter, key=lambda x: (x[1], x[0]), reverse=True):
            if mass >= 1-epsilon or compound.name in zeromass:
                good += 1
            total += 1
            print '{}: {}'.format(compound.translate(lambda x: compound_name.get(x, x)), mass)
        logger.info('Consistent compounds: {}/{}'.format(good, total))

        logger.info('Is consistent? {}'.format(mass_consistency.is_consistent(model, known_inconsistent, zeromass)))

        logger.info('Reaction consistency')
        reaction_iter, compound_iter = mass_consistency.check_reaction_consistency(model, known_inconsistent, zeromass)
        for reaction_id, residual in sorted(reaction_iter, key=lambda x: abs(x[1]), reverse=True):
            if abs(residual) >= epsilon:
                reaction = model.get_reaction(reaction_id).translated_compounds(lambda x: compound_name.get(x, x))
                print '{}\t{}\t{}'.format(reaction_id, residual, reaction)

class RobustnessCommand(Command):
    '''Run robustness analysis on metabolic model

    Given a reaction to maximize and a reaction to vary,
    the robustness analysis will run FBA while fixing the
    reaction to vary at each iteration. The reaction will
    be fixed at the specified number of steps between the
    minimum and maximum flux value specified in the model.'''

    name = 'robustness'
    title = 'Run robustness analysis on a metabolic model'

    def init_parser(self, parser):
        parser.add_argument('--steps', metavar='N', help='Number of flux value steps for varying reaction',
                            type=int, default=10)
        parser.add_argument('--minimum', metavar='V', help='Minumum flux value of varying reacton',
                            type=float)
        parser.add_argument('--maximum', metavar='V', help='Maximum flux value of varying reacton',
                            type=float)
        parser.add_argument('--no-tfba', help='Disable thermodynamic constraints on FBA',
                            action='store_true')
        parser.add_argument('reaction', help='Reaction to maximize')
        parser.add_argument('varying', help='Reaction to vary')

    def __call__(self, model, compounds, **kwargs):
        '''Run flux analysis command'''

        from .lpsolver import cplex
        solver = cplex.Solver()

        def run_fba_fmin(model, reaction):
            fba = fluxanalysis.FluxBalanceProblem(model, solver=solver)
            fba.solve(reaction)
            optimum = fba.get_flux(reaction)
            return fluxanalysis.flux_minimization(model, { reaction: optimum }, solver=solver)

        def run_tfba(model, reaction):
            return fluxanalysis.flux_balance_td(model, reaction, solver=solver)

        if kwargs['no_tfba']:
            run_fba = run_fba_fmin
        else:
            run_fba = run_tfba

        # Load compound information
        compound_name = {}
        for compound in compounds:
            compound_name[compound.id] = compound.name if compound.name is not None else compound.id

        reaction = kwargs['reaction']
        if not model.has_reaction(reaction):
            raise ValueError('Specified reaction is not in model: {}'.format(reaction))

        varying_reaction = kwargs['varying']
        if not model.has_reaction(varying_reaction):
            raise ValueError('Specified reaction is not in model: {}'.format(varying_reaction))

        steps = kwargs['steps']
        if steps <= 0:
            raise ValueError('Invalid number of steps: {}\n'.format(steps))

        # Run FBA on model at different fixed flux values
        flux_min = kwargs['minimum'] if kwargs['minimum'] is not None else model.limits[varying_reaction].lower
        flux_max = kwargs['maximum'] if kwargs['maximum'] is not None else model.limits[varying_reaction].upper

        if flux_min > flux_max:
            raise ValueError('Invalid flux range: {}, {}\n'.format(flux_min, flux_max))

        for i in xrange(steps):
            fixed_flux = flux_min + i*(flux_max - flux_min)/float(steps-1)
            test_model = model.copy()
            test_model.limits[varying_reaction].bounds = fixed_flux, fixed_flux

            try:
                for other_reaction, flux in run_fba(test_model, reaction):
                    print '{}\t{}\t{}'.format(other_reaction, fixed_flux, flux)
            except fluxanalysis.FluxBalanceError:
                pass

class SearchCommand(Command):
    '''Defines the search command'''

    name = 'search'
    title = 'Search the database of reactions or compounds'

    def init_parser(self, parser):
        '''Initialize argument parser'''
        subparsers = parser.add_subparsers(title='Search domain')

        # Compound subcommand
        parser_compound = subparsers.add_parser('compound', help='Search in compounds')
        parser_compound.set_defaults(which='compound')
        parser_compound.add_argument('--id', '-i', dest='id', metavar='id', type=str,
                                        default=[], action='append', help='Compound ID')
        parser_compound.add_argument('--name', '-n', dest='name', metavar='name', type=str,
                                        default=[], action='append', help='Name of compound')

        # Reaction subcommand
        parser_reaction = subparsers.add_parser('reaction', help='Search in reactions')
        parser_reaction.set_defaults(which='reaction')
        parser_reaction.add_argument('--id', '-i', dest='id', metavar='id', type=str,
                                        default=[], action='append', help='Reaction ID')
        parser_reaction.add_argument('--compound', '-c', dest='compound', metavar='compound', type=str,
                                        default=[], action='append', help='Comma-separated list of compound IDs')

    def __call__(self, model, compounds, **kwargs):
        '''Run search command'''

        def parse_compound(s):
            '''Parse a compound specification with optional compartment'''
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
        for compound in compounds:
            compound_name[compound.id] = compound.name if compound.name is not None else compound.id
            compound_synonyms[compound.id] = compound.names

            for n in compound.names:
                n = filter_search_term(n)
                compound_for_name[n] = compound.id

            if compound.formula is not None and not '.' in compound.formula:
                compound_formula[compound.id] = Formula.parse(compound.formula)

        which_command = kwargs['which']
        if which_command == 'compound':
            compound_ids = []
            for name in kwargs['name']:
                search_name = filter_search_term(name)
                if search_name in compound_for_name:
                    compound_ids.append(compound_for_name[search_name])
            for compound_id in kwargs['id']:
                if compound_id in compound_name:
                    compound_ids.append(compound_id)

            # Show results
            for compound_id in compound_ids:
                print 'ID: {}'.format(compound_id)
                if compound_id in compound_name:
                    print 'Name: {}'.format(compound_name[compound_id])
                if compound_id in compound_synonyms:
                    print 'Synonyms: {}'.format(', '.join(compound_synonyms[compound_id]))
                if compound_id in compound_formula:
                    print 'Formula: {}'.format(compound_formula[compound_id])
                print
        elif which_command == 'reaction':
            reaction_ids = []
            for reaction_id in kwargs['id']:
                if reaction_id in model.reactions:
                    reaction_ids.append(reaction_id)
            for compound_list in kwargs['compound']:
                matched_reactions = None
                for compound_spec in compound_list.split(','):
                    compound = parse_compound(compound_spec.strip())
                    if matched_reactions is None:
                        matched_reactions = set(model.get_compound_reactions(compound))
                    else:
                        matched_reactions &= set(model.get_compound_reactions(compound))
                if matched_reactions is not None:
                    reaction_ids.extend(matched_reactions)

            # Show results
            for reaction_id in reaction_ids:
                print 'ID: {}'.format(reaction_id)

                reaction = model.get_reaction(reaction_id)
                print 'Reaction (IDs): {}'.format(model.get_reaction(reaction_id))
                translated_reaction = reaction.translated_compounds(lambda x: compound_name.get(x, x))
                if reaction != translated_reaction:
                    print 'Reaction (names): {}'.format(translated_reaction)
                print


def main(command=None):
    '''Run the command line interface with the given Command'''

    title = 'Metabolic modeling tools'
    if command is not None:
        title = command.title

    parser = argparse.ArgumentParser(description=title)
    parser.add_argument('--database', metavar='file', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Files to use as reaction database')
    parser.add_argument('--compounds', metavar='file', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Files to use as compound database')
    parser.add_argument('--model', metavar='model', nargs=1,
                        type=argparse.FileType('r'),
                        help='File to use as model definition (database subset)')
    parser.add_argument('--limits', metavar='file', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Optional limits on flux of reactions')

    if command is not None:
        # Command explicitly given, only allow that command
        command.init_parser(parser)
        parser.set_defaults(command=command)
    else:
        # Discover all available commands
        commands = {}
        for command_class in Command.__subclasses__():
            command_name = getattr(command_class, 'name', None)
            if command_name is not None and command_name not in commands:
                try:
                    commands[command_name] = command_class.title, command_class()
                except:
                    pass

        # Create parsers for subcommands
        subparsers = parser.add_subparsers(title='Command')
        for name, values in sorted(commands.iteritems()):
            title, command = values
            subparser = subparsers.add_parser(name, help=title)
            subparser.set_defaults(command=command)
            command.init_parser(subparser)

    args = parser.parse_args()
    command = args.command

    # Load reaction database from file
    if len(args.database) == 0:
        raise ValueError('No database provided')

    databases = []
    for database_file in args.database:
        db = DictDatabase()
        for reaction_id, reaction in internal.parse_reaction_file(database_file):
            db.set_reaction(reaction_id, reaction)
        databases.append(db)
    database = ChainedDatabase(*databases)

    if args.model is not None:
        # Set database and model to the database subset
        model = MetabolicModel.load_model(database, internal.parse_model_file(args.model[0]))
    else:
        # Build model from all database reactions
        model = MetabolicModel.load_model(database, database.reactions)

    # Load bounds on exchange reactions
    for limits_table in args.limits:
        for reaction_id, lower, upper in internal.parse_limits_file(limits_table):
            if model.has_reaction(reaction_id):
                if lower is not None:
                    model.limits[reaction_id].lower = lower
                if upper is not None:
                    model.limits[reaction_id].upper = upper

    # Parse compound tables
    def compound_iter():
        for compound_table in args.compounds:
            for compound in modelseed.parse_compound_file(compound_table):
                yield compound

    # Call command
    arg_filter = ('database', 'compounds', 'model', 'limits', 'command')
    kwargs = { key: value for key, value in vars(args).iteritems() if key not in arg_filter }

    command(model=model, compounds=compound_iter(), **kwargs)
