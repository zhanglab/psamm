
'''Utilities for the command line interface'''

import sys
import argparse
import operator
import re

from .fastcore import Fastcore
from .fluxanalysis import flux_balance, flux_minimization
from .formula import Formula, Radical
from .gapfill import gapfind, gapfill
from .massconsistency import MassConsistencyCheck
from .metabolicmodel import DictDatabase
from .reaction import Reaction, Compound
from . import lpsolver, modelseed

class Command(object):
    '''Represents a command in the interface, operating on a model or database

    Subclasses must implement __call__ to handle command execution. Arguments
    will be given as keyword arguments. The keywords database, model and
    compounds will be set.

    In addition, init_parser() can be implemented which will allow the
    command to initialize an instance of ArgumentParser as desired. The
    resulting arguments will be given as keyword arguments to __call__.
    '''

    def __init__(self, title=None):
        self._title = title

    @property
    def title(self):
        '''Command title'''
        return self._title

    def init_parser(self, parser):
        '''Initialize command line parser (argparse.ArgumentParser)'''
        pass

    def __call__(self):
        '''Execute command'''
        pass

class ConsoleCommand(Command):
    '''Start an interactive Python console with the given model loaded'''

    def __init__(self):
        super(ConsoleCommand, self).__init__('Start Python console with metabolic model loaded')

    def __call__(self, **kwargs):
        # Importing readline will in some cases print weird escape
        # characters to stdout. To avoid this we only import readline
        # and related packages at this point when we are certain
        # they are needed.
        from code import InteractiveConsole
        import readline, rlcompleter

        readline.set_completer(rlcompleter.Completer(kwargs).complete)
        readline.parse_and_bind('tab: complete')
        console = InteractiveConsole(kwargs)
        console.interact('Metabolic model has been loaded into "model", "database" and "compounds".')

class FastGapFillCommand(Command):
    def __init__(self):
        super(FastGapFillCommand, self).__init__('Run FastGapFill on a metabolic model')

    def init_parser(self, parser):
        parser.add_argument('--limits', metavar='file', action='append',
                            type=argparse.FileType('r'), default=[],
                            help='Optional limits on flux of reactions')
        parser.add_argument('--penalty', metavar='file', type=argparse.FileType('r'),
                            help='List of penalty scores for database reactions')
        parser.add_argument('reaction', help='Reaction to maximize')

    def __call__(self, database, model, compounds=[], **kwargs):
        '''Run FastGapFill command'''

        if model is None:
            sys.stderr.write('Please specify a model. Use --help for more information.\n')
            sys.exit(-1)

        # Create fastcore object
        solver = lpsolver.CplexSolver(None)
        fastcore = Fastcore(solver)

        # Load compound information
        compound_name = {}
        for compound_table in compounds:
            for compound in modelseed.parse_compound_file(compound_table):
                compound_name[compound.id] = compound.name if compound.name is not None else compound.id

        epsilon = 1e-5
        model_compartments = { None, 'e' }

        # Add exchange and transport reactions to database
        model_complete = model.copy()
        print 'Adding database, exchange and transport reactions...'
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
        print 'Calculating Fastcore induced set on model...'
        core = set(model.reactions)

        induced = fastcore.fastcore(model_complete, core, epsilon, weights=weights)
        print 'Result: |A| = {}, A = {}'.format(len(induced), induced)
        added_reactions = induced - core
        print 'Extended: |E| = {}, E = {}'.format(len(added_reactions), added_reactions)

        # Load bounds on exchange reactions
        for limits_table in kwargs['limits']:
            model.load_reaction_limits(limits_table)

        maximized_reaction = kwargs['reaction']
        print 'Flux balance on induced model maximizing {}...'.format(maximized_reaction)
        model_induced = model.copy()
        for rxnid in induced:
            model_induced.add_reaction(rxnid)
        for rxnid, flux in sorted(flux_balance(model_induced, maximized_reaction)):
            reaction_class = 'Dbase'
            weight = weights.get(rxnid, 1)
            if model.has_reaction(rxnid):
                reaction_class = 'Model'
                weight = 0
            reaction = model_complete.get_reaction(rxnid).translated_compounds(lambda x: compound_name.get(x, x))
            print '{}\t{}\t{}\t{}\t{}'.format(rxnid, reaction_class, weight, flux, reaction)

        print 'Calculating Fastcc consistent subset of induced model...'
        consistent_core = fastcore.fastcc_consistent_subset(model_induced, epsilon)
        print 'Result: |A| = {}, A = {}'.format(len(consistent_core), consistent_core)
        removed_reactions = set(model_induced.reactions) - consistent_core
        print 'Removed: |R| = {}, R = {}'.format(len(removed_reactions), removed_reactions)

class FluxAnalysisCommand(Command):
    def __init__(self):
        super(FluxAnalysisCommand, self).__init__('Run flux analysis on a metabolic model')

    def init_parser(self, parser):
        parser.add_argument('--limits', metavar='limitsfile', action='append',
                            type=argparse.FileType('r'), default=[],
                            help='Optional limits on flux of reactions')
        parser.add_argument('reaction', help='Reaction to maximize')

    def __call__(self, database, model, compounds=[], **kwargs):
        '''Run flux analysis command'''

        if model is None:
            sys.stderr.write('Please specify a model. Use --help for more information.\n')
            sys.exit(-1)

        # Load compound information
        compound_name = {}
        for compound_table in compounds:
            for compound in modelseed.parse_compound_file(compound_table):
                compound_name[compound.id] = compound.name if compound.name is not None else compound.id

        reaction = kwargs['reaction']
        if not model.has_reaction(reaction):
            sys.stderr.write('Specified reaction is not in model: {}\n'.format(reaction))
            sys.exit(-1)

        for limits_table in kwargs['limits']:
            model.load_reaction_limits(limits_table)

        epsilon = 1e-5

        # Run FBA on model
        fba_fluxes = dict(flux_balance(model, reaction))
        optimum = fba_fluxes[reaction]

        # Run flux minimization
        fmin_fluxes = dict(flux_minimization(model, { reaction: optimum }))
        count = 0
        for rxnid, flux in sorted(fmin_fluxes.iteritems()):
            if fba_fluxes[rxnid] - epsilon > flux:
                count += 1
            rx = database.get_reaction(rxnid)
            print '{}\t{}\t{}\t{}'.format(rxnid, fba_fluxes[rxnid], flux,
                                            rx.translated_compounds(lambda x: compound_name.get(x, x)))
        print 'Maximum flux: {}'.format(optimum)
        print 'Minimized reactions: {}'.format(count)

class FormulaBalanceCommand(Command):
    def __init__(self):
        super(FormulaBalanceCommand, self).__init__('Check formula balance on a model or database')

    def __call__(self, database, model, compounds):
        '''Run formula balance command'''

        # Mapping from compound id to formula
        compound_formula = {}
        for compound_table in compounds:
            for compound in modelseed.parse_compound_file(compound_table):
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
                        print 'Error parsing {}: {}'.format(compound.formula, e)
                compound_formula[compound.id] = f

        # Create a set of known mass-inconsistent reactions
        exchange = set()
        for rxnid in database.reactions:
            rx = database.get_reaction(rxnid)
            if len(rx.left) == 0 or len(rx.right) == 0:
                exchange.add(rxnid)

        def multiply_formula(compound_list):
            for compound, count in compound_list:
                yield count * compound_formula.get(compound.name, Formula())

        for reaction in database.reactions:
            if reaction not in exchange:
                rx = database.get_reaction(reaction)
                left_form = reduce(operator.or_, multiply_formula(rx.left), Formula())
                right_form = reduce(operator.or_, multiply_formula(rx.right), Formula())

                if right_form != left_form:
                    right_missing, left_missing = Formula.balance(right_form, left_form)
                    print '{}\t{}\t{}\t{}\t{}'.format(reaction, left_form, right_form, left_missing, right_missing)

class GapFillCommand(Command):
    def __init__(self):
        super(GapFillCommand, self).__init__('Run GapFind and GapFill on a metabolic model')

    def __call__(self, database, model, compounds=[]):
        '''Run GapFill command'''

        if model is None:
            sys.stderr.write('Please specify a model. Use --help for more information.\n')
            sys.exit(-1)

        # Load compound information
        compound_name = {}
        for compound_table in compounds:
            for compound in modelseed.parse_compound_file(compound_table):
                compound_name[compound.id] = compound.name if compound.name is not None else compound.id

        model_compartments = { None, 'e' }

        # Run GapFind on model
        print 'Searching for blocked compounds...'
        blocked = set(compound for compound in gapfind(model) if compound.compartment is not 'e')
        if len(blocked) > 0:
            print 'Blocked:'
            for compound in blocked:
                print compound.translate(lambda x: compound_name.get(x, x))

        if len(blocked) > 0:
            # Add exchange and transport reactions to database
            model_complete = model.copy()
            print 'Adding database, exchange and transport reactions...'
            model_complete.add_all_database_reactions(model_compartments)
            model_complete.add_all_exchange_reactions()
            model_complete.add_all_transport_reactions()

            print 'Searching for reactions to fill gaps...'
            added_reactions, reversed_reactions = gapfill(model_complete, model.reactions, blocked)

            for rxnid in added_reactions:
                rx = model_complete.get_reaction(rxnid).translated_compounds(lambda x: compound_name.get(x, x))
                print '{}\t{}\t{}'.format(rxnid, 'Add', rx)

            for rxnid in reversed_reactions:
                rx = model_complete.get_reaction(rxnid).translated_compounds(lambda x: compound_name.get(x, x))
                print '{}\t{}\t{}'.format(rxnid, 'Reverse', rx)
        else:
            print 'No blocked compounds found'

class MassConsistencyCommand(Command):
    def __init__(self):
        super(MassConsistencyCommand, self).__init__('Run mass consistency check on a database')

    def init_parser(self, parser):
        parser.add_argument('--exclude', metavar='reaction', action='append',
                            type=str, default=[], help='Exclude reaction from mass consistency')

    def __call__(self, database, model, compounds=[], **kwargs):
        # Load compound information
        compound_name = {}
        for compound_table in compounds:
            for compound in modelseed.parse_compound_file(compound_table):
                compound_name[compound.id] = compound.name if compound.name is not None else compound.id

        # Create a set of known mass-inconsistent reactions
        exchange = set()
        for reaction_id in database.reactions:
            rx = database.get_reaction(reaction_id)
            if len(rx.left) == 0 or len(rx.right) == 0:
                exchange.add(reaction_id)

        # Other reactions to exclude from consistency check
        exclude = set(kwargs['exclude'])

        # Create set of compounds allowed to have mass zero
        zeromass = set()
        zeromass.add('cpd11632') # Photon
        zeromass.add('cpd12713') # Electron

        mass_consistency = MassConsistencyCheck()
        known_inconsistent = exclude | exchange

        print 'Mass consistency on database...'
        epsilon = 1e-5
        compound_iter = mass_consistency.check_compound_consistency(database, known_inconsistent, zeromass)

        print 'Compound consistency...'
        good = 0
        total = 0
        for compound, mass in sorted(compound_iter, key=lambda x: (x[1], x[0]), reverse=True):
            if mass >= 1-epsilon or compound.name in zeromass:
                good += 1
            total += 1
            print '{}: {}'.format(compound.translate(lambda x: compound_name.get(x, x)), mass)
        print 'Consistent compounds: {}/{}'.format(good, total)

        print 'Is consistent? {}'.format(mass_consistency.is_consistent(database, known_inconsistent, zeromass))

        print 'Reaction consistency...'
        reaction_iter, compound_iter = mass_consistency.check_reaction_consistency(database, known_inconsistent, zeromass)
        for reaction_id, residual in sorted(reaction_iter, key=lambda x: abs(x[1]), reverse=True):
            if abs(residual) >= epsilon:
                reaction = database.get_reaction(reaction_id).translated_compounds(lambda x: compound_name.get(x, x))
                print '{}\t{}\t{}'.format(reaction_id, residual, reaction)

class RobustnessCommand(Command):
    '''Run robustness analysis on metabolic model

    Given a reaction to maximize and a reaction to vary,
    the robustness analysis will run FBA while fixing the
    reaction to vary at each iteration. The reaction will
    be fixed at the specified number of steps between the
    minimum and maximum flux value specified in the model.'''

    def __init__(self):
        super(RobustnessCommand, self).__init__('Run robustness analysis on a metabolic model')

    def init_parser(self, parser):
        parser.add_argument('--limits', metavar='limitsfile', action='append',
                            type=argparse.FileType('r'), default=[],
                            help='Optional limits on flux of reactions')
        parser.add_argument('--steps', metavar='N', help='Number of flux value steps for varying reaction',
                            type=int, default=10)
        parser.add_argument('--minimum', metavar='V', help='Minumum flux value of varying reacton',
                            type=float)
        parser.add_argument('--maximum', metavar='V', help='Maximum flux value of varying reacton',
                            type=float)
        parser.add_argument('reaction', help='Reaction to maximize')
        parser.add_argument('varying', help='Reaction to vary')

    def __call__(self, database, model, compounds=[], **kwargs):
        '''Run flux analysis command'''

        if model is None:
            sys.stderr.write('Please specify a model. Use --help for more information.\n')
            sys.exit(-1)

        # Load compound information
        compound_name = {}
        for compound_table in compounds:
            for compound in modelseed.parse_compound_file(compound_table):
                compound_name[compound.id] = compound.name if compound.name is not None else compound.id

        reaction = kwargs['reaction']
        if not model.has_reaction(reaction):
            sys.stderr.write('Specified reaction is not in model: {}\n'.format(reaction))
            sys.exit(-1)

        varying_reaction = kwargs['varying']
        if not model.has_reaction(varying_reaction):
            sys.stderr.write('Specified reaction is not in model: {}\n'.format(varying_reaction))
            sys.exit(-1)

        steps = kwargs['steps']
        if steps <= 0:
            sys.stderr.write('Invalid number of steps: {}\n'.format(steps))
            sys.exit(-1)

        # Run FBA on model at different fixed flux values
        flux_min = kwargs.get('minimum', model.limits[varying_reaction].lower)
        flux_max = kwargs.get('maximum', model.limits[varying_reaction].upper)

        if flux_min > flux_max:
            sys.stderr.write('Invalid flux range: {}, {}\n'.format(flux_min, flux_max))
            sys.exit(-1)

        # Load reaction limits
        for limits_table in kwargs['limits']:
            model.load_reaction_limits(limits_table)

        for i in xrange(steps):
            fixed_flux = flux_min + i*(flux_max - flux_min)/float(steps-1)
            test_model = model.copy()
            test_model.limits[varying_reaction].bounds = fixed_flux, fixed_flux

            try:
                fba_fluxes = dict(flux_balance(test_model, reaction))
                optimum = fba_fluxes[reaction]
                fmin_fluxes = dict(flux_minimization(test_model, { reaction: optimum }))
                for other_reaction, flux in fmin_fluxes.iteritems():
                    print '{}\t{}\t{}'.format(other_reaction, fixed_flux, flux)
            except Exception:
                pass

class SearchCommand(Command):
    '''Defines the search command'''

    def __init__(self):
        super(SearchCommand, self).__init__('Search the database of reactions or compounds')

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

    def __call__(self, database, model, compounds=[], **kwargs):
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
        for compound_table in compounds:
            for compound in modelseed.parse_compound_file(compound_table):
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
                if reaction_id in database.reactions:
                    reaction_ids.append(reaction_id)
            for compound_list in kwargs['compound']:
                matched_reactions = None
                for compound_spec in compound_list.split(','):
                    compound = parse_compound(compound_spec.strip())
                    if matched_reactions is None:
                        matched_reactions = set(database.get_compound_reactions(compound))
                    else:
                        matched_reactions &= database.compound_reactions.get(compound, set())
                if matched_reactions is not None:
                    reaction_ids.extend(matched_reactions)

            # Show results
            for reaction_id in reaction_ids:
                print 'ID: {}'.format(reaction_id)

                reaction = database.get_reaction(reaction_id)
                print 'Reaction (IDs): {}'.format(database.get_reaction(reaction_id))
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

    if command is not None:
        # Command explicitly given, only allow that command
        command.init_parser(parser)
        parser.set_defaults(command=command)
    else:
        # Allow all commands as subcommands
        subcommands = {
            'console': ConsoleCommand(),
            'fastgapfill': FastGapFillCommand(),
            'fba': FluxAnalysisCommand(),
            'formulacheck': FormulaBalanceCommand(),
            'gapfill': GapFillCommand(),
            'masscheck': MassConsistencyCommand(),
            'robustness': RobustnessCommand(),
            'search': SearchCommand()
        }

        subparsers = parser.add_subparsers(title='Command')
        for name, command in sorted(subcommands.iteritems()):
            subparser = subparsers.add_parser(name, help=command.title)
            subparser.set_defaults(command=command)
            command.init_parser(subparser)

    args = parser.parse_args()
    command = args.command

    # Load reaction database from file
    database = DictDatabase.load_from_files(*args.database)
    model = None
    if args.model is not None:
        # Set database and model to the database subset
        database = database.load_model_from_file(args.model[0])
        model = database

    compounds = args.compounds

    # Call command
    arg_filter = ('database', 'compounds', 'model', 'command')
    kwargs = { key: value for key, value in vars(args).iteritems() if key not in arg_filter }
    command(database=database, model=model, compounds=compounds, **kwargs)
