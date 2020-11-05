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
# Copyright 2015-2017  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>

"""Entry points and functionality for importing models.

This module contains entry points for the psamm-import programs and
functionality to assist in importing external file formats, converting
them into proper native models and writing those models to files.
"""

from __future__ import print_function, unicode_literals

import sys
import os
import argparse
import logging
import math
import re
import json
import codecs
from itertools import product
from collections import OrderedDict, Counter
from decimal import Decimal

import yaml
import pkg_resources
from six import iteritems, itervalues, text_type
from six.moves.urllib.request import urlopen
from six.moves.urllib.parse import quote as url_quote

from .datasource.native import ModelWriter
from .datasource.reaction import (parse_reaction,
                                  ParseError as ReactionParseError)
from .datasource.entry import DictCompartmentEntry
from .datasource import sbml
from .expression import boolean
from . import formula

from .util import mkdir_p

# Threshold for putting reactions into subsystem files
_MAX_REACTION_COUNT = 3

# Threshold for converting reactions into dictionary representation.
_MAX_REACTION_LENGTH = 10

logger = logging.getLogger(__name__)


class ImportError(Exception):
    """Exception used to signal a general import error."""


class ModelLoadError(ImportError):
    """Exception used to signal an error loading the model files."""


class ParseError(ImportError):
    """Exception used to signal an error parsing the model files."""


class Importer(object):
    """Base importer class."""

    def _try_parse_formula(self, compound_id, s):
        """Try to parse the given compound formula string.

        Logs a warning if the formula could not be parsed.
        """
        s = s.strip()
        if s == '':
            return None

        try:
            # Do not return the parsed formula. For now it is better to keep
            # the original formula string unchanged in all cases.
            formula.Formula.parse(s)
        except formula.ParseError:
            logger.warning('Unable to parse compound formula {}: {}'.format(
                compound_id, s))

        return s

    def _try_parse_reaction(self, reaction_id, s,
                            parser=parse_reaction, **kwargs):
        """Try to parse the given reaction equation string.

        Returns the parsed Reaction object, or raises an error if the reaction
        could not be parsed.
        """
        try:
            return parser(s, **kwargs)
        except ReactionParseError as e:
            if e.indicator is not None:
                logger.error('{}\n{}\n{}'.format(
                    str(e), s, e.indicator))
            raise ParseError('Unable to parse reaction {}: {}'.format(
                reaction_id, s))

    def _try_parse_gene_association(self, reaction_id, s):
        """Try to parse the given gene association rule.

        Logs a warning if the association rule could not be parsed and returns
        the original string. Otherwise, returns the boolean.Expression object.
        """
        s = s.strip()
        if s == '':
            return None

        try:
            return boolean.Expression(s)
        except boolean.ParseError as e:
            msg = 'Failed to parse gene association for {}: {}'.format(
                reaction_id, text_type(e))
            if e.indicator is not None:
                msg += '\n{}\n{}'.format(s, e.indicator)
            logger.warning(msg)

        return s


# Define custom dict representers for YAML
# This allows reading/writing Python OrderedDicts in the correct order.
# See: https://stackoverflow.com/questions/5121931/in-python-how-can-you-load-yaml-mappings-as-ordereddicts  # noqa
def _dict_representer(dumper, data):
    return dumper.represent_dict(iteritems(data))


def _dict_constructor(loader, node):
    return OrderedDict(loader.construct_pairs(node))


def _set_representer(dumper, data):
    return dumper.represent_list(iter(data))


def _boolean_expression_representer(dumper, data):
    return dumper.represent_unicode(text_type(data))


def _decimal_representer(dumper, data):
    # Code from float_representer in PyYAML.
    if data % 1 == 0:
        return dumper.represent_int(int(data))
    elif math.isnan(data):
        value = '.nan'
    elif data == float('inf'):
        value = '.inf'
    elif data == float('-inf'):
        value = '-.inf'
    else:
        value = text_type(data).lower()
        if '.' not in value and 'e' in value:
            value = value.replace('e', '.0e', 1)
    return dumper.represent_scalar('tag:yaml.org,2002:float', value)


def get_default_compartment(model):
    """Return what the default compartment should be set to.

    If some compounds have no compartment, unique compartment
    name is returned to avoid collisions.
    """
    default_compartment = 'c'
    default_key = set()
    for reaction in model.reactions:
        equation = reaction.equation
        if equation is None:
            continue

        for compound, _ in equation.compounds:
            default_key.add(compound.compartment)

    if None in default_key and default_compartment in default_key:
        suffix = 1
        while True:
            new_key = '{}_{}'.format(default_compartment, suffix)
            if new_key not in default_key:
                default_compartment = new_key
                break
            suffix += 1

    if None in default_key:
        logger.warning(
            'Compound(s) found without compartment, default'
            ' compartment is set to {}.'.format(default_compartment))
    return default_compartment


def detect_best_flux_limit(model):
    """Detect the best default flux limit to use for model output.

    The default flux limit does not change the model but selecting a good
    value reduced the amount of output produced and reduces clutter in the
    output files.
    """
    flux_limit_count = Counter()

    for reaction in model.reactions:
        if reaction.id not in model.limits:
            continue

        equation = reaction.properties['equation']
        if equation is None:
            continue

        _, lower, upper = model.limits[reaction.id]
        if upper is not None and upper > 0 and equation.direction.forward:
            flux_limit_count[upper] += 1
        if lower is not None and -lower > 0 and equation.direction.reverse:
            flux_limit_count[-lower] += 1

    if len(flux_limit_count) == 0:
        return None

    best_flux_limit, _ = flux_limit_count.most_common(1)[0]
    return best_flux_limit


def reactions_to_files(model, dest, writer, split_subsystem):
    """Turn the reaction subsystems into their own files.

    If a subsystem has a number of reactions over the threshold, it gets its
    own YAML file. All other reactions, those that don't have a subsystem or
    are in a subsystem that falls below the threshold, get added to a common
    reaction file.

    Args:
        model: :class:`psamm_import.model.MetabolicModel`.
        dest: output path for model files.
        writer: :class:`psamm.datasource.native.ModelWriter`.
        split_subsystem: Divide reactions into multiple files by subsystem.
    """
    def safe_file_name(origin_name):
        safe_name = re.sub(
            r'\W+', '_', origin_name, flags=re.UNICODE)
        safe_name = re.sub(
            r'_+', '_', safe_name.lower(), flags=re.UNICODE)
        safe_name = safe_name.strip('_')
        return safe_name

    common_reactions = []
    reaction_files = []
    if not split_subsystem:
        common_reactions = sorted(model.reactions, key=lambda r: r.id)
        if len(common_reactions) > 0:
            reaction_file = 'reactions.yaml'
            with open(os.path.join(dest, reaction_file), 'w') as f:
                writer.write_reactions(f, common_reactions)
            reaction_files.append(reaction_file)
    else:
        subsystems = {}
        for reaction in sorted(model.reactions, key=lambda r: r.id):
            if 'subsystem' in reaction.properties:
                subsystem_file = safe_file_name(
                    reaction.properties['subsystem'])
                subsystems.setdefault(subsystem_file, []).append(reaction)
            else:
                common_reactions.append(reaction)

        subsystem_folder = 'reactions'
        sub_existance = False
        for subsystem_file, reactions in iteritems(subsystems):
            if len(reactions) < _MAX_REACTION_COUNT:
                for reaction in reactions:
                    common_reactions.append(reaction)
            else:
                if len(reactions) > 0:
                    mkdir_p(os.path.join(dest, subsystem_folder))
                    subsystem_file = os.path.join(
                        subsystem_folder, '{}.yaml'.format(subsystem_file))

                    with open(os.path.join(dest, subsystem_file), 'w') as f:
                        writer.write_reactions(f, reactions)
                    reaction_files.append(subsystem_file)
                    sub_existance = True

        reaction_files.sort()
        if sub_existance:
            reaction_file = os.path.join(
                subsystem_folder, 'other_reactions.yaml')
        else:
            reaction_file = 'reactions.yaml'
        if len(common_reactions) > 0:
            with open(os.path.join(dest, reaction_file), 'w') as f:
                writer.write_reactions(f, common_reactions)
            reaction_files.append(reaction_file)

    return reaction_files


def _get_output_limit(limit, default_limit):
    """Return limit to output given default limit."""
    return limit if limit != default_limit else None


def _generate_limit_items(lower, upper):
    """Yield key, value pairs for limits dictionary.

    Yield pairs of key, value where key is ``lower``, ``upper`` or ``fixed``.
    A key, value pair is emitted if the bounds are not None.
    """
    # Use value + 0 to convert any -0.0 to 0.0 which looks better.
    if lower is not None and upper is not None and lower == upper:
        yield 'fixed', upper + 0
    else:
        if lower is not None:
            yield 'lower', lower + 0
        if upper is not None:
            yield 'upper', upper + 0


def model_exchange(model):
    """Return exchange definition as YAML dict."""
    # Determine the default flux limits. If the value is already at the
    # default it does not need to be included in the output.
    lower_default, upper_default = None, None
    if model.default_flux_limit is not None:
        lower_default = -model.default_flux_limit
        upper_default = model.default_flux_limit

    compounds = []
    for compound, reaction_id, lower, upper in sorted(
            itervalues(model.exchange)):
        d = OrderedDict([('id', compound.name)])
        if reaction_id is not None:
            d['reaction'] = reaction_id

        lower = _get_output_limit(lower, lower_default)
        upper = _get_output_limit(upper, upper_default)
        d.update(_generate_limit_items(lower, upper))

        compounds.append(d)

    return OrderedDict([('compounds', compounds)])


def model_reaction_limits(model):
    """Yield model reaction limits as YAML dicts."""
    for reaction in sorted(model.reactions, key=lambda r: r.id):
        equation = reaction.properties.get('equation')
        if equation is None:
            continue

        # Determine the default flux limits. If the value is already at the
        # default it does not need to be included in the output.
        lower_default, upper_default = None, None
        if model.default_flux_limit is not None:
            if equation.direction.reverse:
                lower_default = -model.default_flux_limit
            else:
                lower_default = 0.0

            if equation.direction.forward:
                upper_default = model.default_flux_limit
            else:
                upper_default = 0.0

        lower_flux, upper_flux = None, None
        if reaction.id in model.limits:
            _, lower, upper = model.limits[reaction.id]
            lower_flux = _get_output_limit(lower, lower_default)
            upper_flux = _get_output_limit(upper, upper_default)

        if lower_flux is not None or upper_flux is not None:
            d = OrderedDict([('reaction', reaction.id)])
            d.update(_generate_limit_items(lower_flux, upper_flux))

            yield d


def infer_compartment_entries(model):
    """Infer compartment entries for model based on reaction compounds."""
    compartment_ids = set()
    for reaction in model.reactions:
        equation = reaction.equation
        if equation is None:
            continue

        for compound, _ in equation.compounds:
            compartment = compound.compartment
            if compartment is None:
                compartment = model.default_compartment

            if compartment is not None:
                compartment_ids.add(compartment)

    for compartment in compartment_ids:
        if compartment in model.compartments:
            continue

        entry = DictCompartmentEntry(dict(id=compartment))
        model.compartments.add_entry(entry)


def infer_compartment_adjacency(model):
    """Infer compartment adjacency for model based on reactions."""
    def reaction_compartments(seq):
        for compound, _ in seq:
            compartment = compound.compartment
            if compartment is None:
                compartment = model.default_compartment

            if compartment is not None:
                yield compartment

    for reaction in model.reactions:
        equation = reaction.equation
        if equation is None:
            continue

        left = reaction_compartments(equation.left)
        right = reaction_compartments(equation.right)
        for c1, c2 in product(left, right):
            if c1 == c2:
                continue
            model.compartment_boundaries.add((c1, c2))
            model.compartment_boundaries.add((c2, c1))


def count_genes(model):
    """Count the number of distinct genes in model reactions."""
    genes = set()
    for reaction in model.reactions:
        if reaction.genes is None:
            continue

        if isinstance(reaction.genes, boolean.Expression):
            genes.update(v.symbol for v in reaction.genes.variables)
        else:
            gene_exp = boolean.Expression(reaction.genes)
            genes.update(v.symbol for v in gene_exp.variables)
            # genes.update(reaction.genes)

    return len(genes)


def write_yaml_model(model, dest='.', convert_exchange=True,
                     split_subsystem=True):
    """Write the given NativeModel to YAML files in dest folder.

    The parameter ``convert_exchange`` indicates whether the exchange reactions
    should be converted automatically to an exchange file.
    """
    yaml.SafeDumper.add_representer(OrderedDict, _dict_representer)
    yaml.SafeDumper.add_representer(set, _set_representer)
    yaml.SafeDumper.add_representer(frozenset, _set_representer)
    yaml.SafeDumper.add_representer(
        boolean.Expression, _boolean_expression_representer)
    yaml.SafeDumper.add_representer(Decimal, _decimal_representer)

    yaml.SafeDumper.ignore_aliases = lambda *args: True

    yaml.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
                         _dict_constructor)

    yaml_args = {'default_flow_style': False,
                 'encoding': 'utf-8',
                 'allow_unicode': True,
                 'width': 79}

    # The ModelWriter from PSAMM is not yet able to write the full model but
    # only reactions and compounds.
    writer = ModelWriter()

    with open(os.path.join(dest, 'compounds.yaml'), 'w+') as f:
        writer.write_compounds(f, sorted(model.compounds, key=lambda c: c.id))

    if model.default_flux_limit is None:
        model.default_flux_limit = detect_best_flux_limit(model)

    if model.extracellular_compartment is None:
        model.extracellular_compartment = (
            sbml.detect_extracellular_compartment(model))

    if model.default_compartment is None:
        model.default_compartment = get_default_compartment(model)

    if model.default_flux_limit is not None:
        logger.info('Using default flux limit of {}'.format(
            model.default_flux_limit))

    if convert_exchange:
        logger.info('Converting exchange reactions to exchange file')
        sbml.convert_exchange_to_compounds(model)

    if len(model.compartments) == 0:
        infer_compartment_entries(model)
        logger.info('Inferred {} compartments: {}'.format(
            len(model.compartments),
            ', '.join(c.id for c in model.compartments)))

    if len(model.compartments) != 0 and len(model.compartment_boundaries) == 0:
        infer_compartment_adjacency(model)

    reaction_files = reactions_to_files(model, dest, writer, split_subsystem)

    if len(model.exchange) > 0:
        with open(os.path.join(dest, 'exchange.yaml'), 'w+') as f:
            yaml.safe_dump(model_exchange(model), f, **yaml_args)

    reaction_limits = list(model_reaction_limits(model))
    if len(reaction_limits) > 0:
        with open(os.path.join(dest, 'limits.yaml'), 'w+') as f:
            yaml.safe_dump(reaction_limits, f, **yaml_args)

    model_d = OrderedDict()
    if model.name is not None:
        model_d['name'] = model.name

    if model.biomass_reaction is not None:
        model_d['biomass'] = model.biomass_reaction
    if model.default_flux_limit is not None:
        model_d['default_flux_limit'] = model.default_flux_limit
    if model.extracellular_compartment != 'e':
        model_d['extracellular'] = model.extracellular_compartment
    if model.default_compartment != 'c':
        model_d['default_compartment'] = model.default_compartment

    if len(model.compartments) > 0:
        adjacency = {}
        for c1, c2 in model.compartment_boundaries:
            adjacency.setdefault(c1, set()).add(c2)
            adjacency.setdefault(c2, set()).add(c1)

        compartment_list = []
        for compartment in sorted(model.compartments, key=lambda c: c.id):
            adjacent = adjacency.get(compartment.id)
            if adjacent is not None and len(adjacent) == 1:
                adjacent = next(iter(adjacent))
            compartment_list.append(writer.convert_compartment_entry(
                compartment, adjacent))

        model_d['compartments'] = compartment_list

    model_d['compounds'] = [{'include': 'compounds.yaml'}]
    model_d['reactions'] = []
    for reaction_file in reaction_files:
        model_d['reactions'].append({'include': reaction_file})

    if len(model.exchange) > 0:
        model_d['exchange'] = [{'include': 'exchange.yaml'}]

    if len(reaction_limits) > 0:
        model_d['limits'] = [{'include': 'limits.yaml'}]

    with open(os.path.join(dest, 'model.yaml'), 'w+') as f:
        yaml.safe_dump(model_d, f, **yaml_args)


def main(importer_class=None, args=None):
    """Entry point for import program.

    If the ``args`` are provided, these should be a list of strings that will
    be used instead of ``sys.argv[1:]``. This is mostly useful for testing.
    """
    parser = argparse.ArgumentParser(
        description='Import from external model formats')
    parser.add_argument('--source', metavar='path', default='.',
                        help='Source directory or file')
    parser.add_argument('--dest', metavar='path', default='.',
                        help='Destination directory (default is ".")')
    parser.add_argument('--no-exchange', action='store_true',
                        help=('Disable importing exchange reactions as'
                              ' exchange compound file.'))
    parser.add_argument('--split-subsystem', action='store_true',
                        help='Enable splitting reaction files by subsystem')
    parser.add_argument('--merge-compounds', action='store_true',
                        help=('Merge identical compounds occuring in various'
                              ' compartments.'))
    parser.add_argument('--force', action='store_true',
                        help='Enable overwriting model files')
    parser.add_argument('--model-name', type=str,
                        help=('specify model name if multiple models '
                              'are stored in the .mat file'))

    if importer_class is None:
        parser.add_argument(
            'format', help='Format to import ("list" to see all)')

    args = parser.parse_args(args)

    # Set up logging for the command line interface
    if 'PSAMM_DEBUG' in os.environ:
        level = getattr(logging, os.environ['PSAMM_DEBUG'].upper(), None)
        if level is not None:
            logging.basicConfig(level=level)
    else:
        logging.basicConfig(
            level=logging.INFO, format='%(levelname)s: %(message)s')

    if importer_class is None:
        # Discover all available model importers
        importers = {}
        for importer_entry in pkg_resources.iter_entry_points(
                'psamm.importer'):
            canonical = importer_entry.name.lower()
            if canonical not in importers:
                importers[canonical] = importer_entry
            else:
                logger.warning('Importer {} was found more than once!'.format(
                    importer_entry.name))

        # Print list of importers
        if args.format in ('list', 'help'):
            if len(importers) == 0:
                logger.error('No importers found!')
            else:
                importer_classes = []
                for name, entry in iteritems(importers):
                    importer_class = entry.load()
                    title = getattr(importer_class, 'title', None)
                    generic = getattr(importer_class, 'generic', False)
                    if title is not None:
                        importer_classes.append(
                            (title, generic, name, importer_class))

                print('Generic importers:')
                for title, _, name, importer_class in sorted(
                        c for c in importer_classes if c[1]):
                    print('{:<12}  {}'.format(name, title))

                print()
                print('Model-specific importers:')
                for title, _, name, importer_class in sorted(
                        c for c in importer_classes if not c[1]):
                    print('{:<12}  {}'.format(name, title))
            sys.exit(0)

        importer_name = args.format.lower()
        if importer_name not in importers:
            logger.error('Importer {} not found!'.format(importer_name))
            logger.info('Use "list" to see available importers.')
            sys.exit(-1)

        importer_class = importers[importer_name].load()

    importer = importer_class()

    try:
        if importer.name == 'matlab' and args.model_name is not None:
            model = importer.import_model(args.source, args.model_name)
        else:
            model = importer.import_model(args.source)
    except ModelLoadError as e:
        logger.error('Failed to load model!', exc_info=True)
        importer.help()
        parser.error(text_type(e))
    except ParseError as e:
        logger.error('Failed to parse model!', exc_info=True)
        logger.error(text_type(e))
        sys.exit(-1)
    except Exception:
        importer.help()

    if args.merge_compounds:
        compounds_before = len(model.compounds)
        sbml.merge_equivalent_compounds(model)
        if len(model.compounds) < compounds_before:
            logger.info(
                'Merged {} compound entries into {} entries by'
                ' removing duplicates in various compartments'.format(
                    compounds_before, len(model.compounds)))

    print('Model: {}'.format(model.name))
    print('- Biomass reaction: {}'.format(model.biomass_reaction))
    print('- Compartments: {}'.format(len(model.compartments)))
    print('- Compounds: {}'.format(len(model.compounds)))
    print('- Reactions: {}'.format(len(model.reactions)))
    print('- Genes: {}'.format(count_genes(model)))

    # Check if dest directory is empty. If we get an error assume that the
    # directory does not exist.
    dest_is_empty = False
    try:
        dest_is_empty = len(os.listdir(args.dest)) == 0
    except OSError:
        dest_is_empty = True

    if not dest_is_empty:
        if not args.force:
            logger.error('Destination directory is not empty. Use --force'
                         ' option to proceed anyway, overwriting any existing'
                         ' files in {}'.format(args.dest))
            return 1
        else:
            logger.warning('Destination directory is not empty, overwriting'
                           ' existing files in {}'.format(args.dest))

    # Create destination directory if not exists
    dest = args.dest
    mkdir_p(dest)

    convert_exchange = not args.no_exchange
    write_yaml_model(model, dest, convert_exchange=convert_exchange,
                     split_subsystem=args.split_subsystem)


def main_bigg(args=None, urlopen=urlopen):
    """Entry point for BiGG import program.

    If the ``args`` are provided, these should be a list of strings that will
    be used instead of ``sys.argv[1:]``. This is mostly useful for testing.
    """
    parser = argparse.ArgumentParser(
        description='Import from BiGG database')
    parser.add_argument('--dest', metavar='path', default='.',
                        help='Destination directory (default is ".")')
    parser.add_argument('--no-exchange', action='store_true',
                        help=('Disable importing exchange reactions as'
                              ' exchange compound file.'))
    parser.add_argument('--split-subsystem', action='store_true',
                        help='Enable splitting reaction files by subsystem')
    parser.add_argument('--merge-compounds', action='store_true',
                        help=('Merge identical compounds occuring in various'
                              ' compartments.'))
    parser.add_argument('--force', action='store_true',
                        help='Enable overwriting model files')
    parser.add_argument('id', help='BiGG model to import ("list" to see all)')

    args = parser.parse_args(args)

    # Set up logging for the command line interface
    if 'PSAMM_DEBUG' in os.environ:
        level = getattr(logging, os.environ['PSAMM_DEBUG'].upper(), None)
        if level is not None:
            logging.basicConfig(level=level)
    else:
        logging.basicConfig(
            level=logging.INFO, format='%(levelname)s: %(message)s')

    # Print list of available models
    if args.id == 'list':
        print('Available models:')
        f = urlopen('http://bigg.ucsd.edu/api/v2/models')
        doc = json.loads(f.read().decode('utf-8'))
        results = doc['results']
        id_width = min(max(len(result['bigg_id']) for result in results), 16)
        for result in sorted(results, key=lambda x: x.get('organism')):
            print('{} {}'.format(
                result.get('bigg_id').ljust(id_width), result.get('organism')))
        return 0

    importer_entry = None
    try:
        importer_entry = next(
            pkg_resources.iter_entry_points('psamm.importer', 'JSON'))
    except StopIteration:
        logger.error('Failed to locate the COBRA JSON model importer!')
        sys.exit(-1)

    importer_class = importer_entry.load()
    importer = importer_class()

    try:
        f = urlopen(
            'http://bigg.ucsd.edu/api/v2/models/{}/download'.format(
                url_quote(args.id)))
        model = importer.import_model(codecs.getreader('utf-8')(f))
    except ModelLoadError as e:
        logger.error('Failed to load model!', exc_info=True)
        importer.help()
        parser.error(text_type(e))
    except ParseError as e:
        logger.error('Failed to parse model!', exc_info=True)
        logger.error(text_type(e))
        sys.exit(-1)

    if args.merge_compounds:
        compounds_before = len(model.compounds)
        sbml.merge_equivalent_compounds(model)
        if len(model.compounds) < compounds_before:
            logger.info(
                'Merged {} compound entries into {} entries by'
                ' removing duplicates in various compartments'.format(
                    compounds_before, len(model.compounds)))

    print('Model: {}'.format(model.name))
    print('- Biomass reaction: {}'.format(model.biomass_reaction))
    print('- Compartments: {}'.format(len(model.compartments)))
    print('- Compounds: {}'.format(len(model.compounds)))
    print('- Reactions: {}'.format(len(model.reactions)))
    print('- Genes: {}'.format(count_genes(model)))

    # Check if dest directory is empty. If we get an error assume that the
    # directory does not exist.
    dest_is_empty = False
    try:
        dest_is_empty = len(os.listdir(args.dest)) == 0
    except OSError:
        dest_is_empty = True

    if not dest_is_empty:
        if not args.force:
            logger.error('Destination directory is not empty. Use --force'
                         ' option to proceed anyway, overwriting any existing'
                         ' files in {}'.format(args.dest))
            return 1
        else:
            logger.warning('Destination directory is not empty, overwriting'
                           ' existing files in {}'.format(args.dest))

    # Create destination directory if not exists
    dest = args.dest
    mkdir_p(dest)

    convert_exchange = not args.no_exchange
    write_yaml_model(model, dest, convert_exchange=convert_exchange,
                     split_subsystem=args.split_subsystem)
