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
# Copyright 2014-2017  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>
# Copyright 2020-2020  Elysha Sameth <esameth1@my.uri.edu>

"""Module for reading and writing native formats.

These formats are either table-based or YAML-based. Table-based formats
are space-separated and empty lines are ignored. Comments starting with
pound (#). YAML-based formats are structured data following the YAML
specification.
"""

from __future__ import absolute_import, unicode_literals

import os
import logging
import re
import csv
import math
from collections import OrderedDict

import yaml
from six import string_types, text_type, iteritems, itervalues, PY3
from decimal import Decimal

from ..reaction import Reaction, Compound, Direction
from ..expression import boolean
from ..formula import Formula
from ..metabolicmodel import MetabolicModel
from ..database import DictDatabase
from .context import FilePathContext, FileMark
from .entry import (DictCompoundEntry as CompoundEntry,
                    DictReactionEntry as ReactionEntry,
                    DictCompartmentEntry as CompartmentEntry)
from .reaction import ReactionParser, convert_to_unicode
from . import modelseed
from .. import util

# Module-level logging
logger = logging.getLogger(__name__)

_HAS_YAML_LIBRARY = None

_REACTION_PARSER = ReactionParser()


class ParseError(Exception):
    """Exception used to signal errors while parsing"""


def float_constructor(loader, node):
    """Construct Decimal from YAML float encoding."""
    s = loader.construct_scalar(node)
    if s == '.inf':
        return Decimal('Infinity')
    elif s == '-.inf':
        return -Decimal('Infinity')
    elif s == '.nan':
        return Decimal('NaN')
    return Decimal(s)


def yaml_load(stream):
    """Load YAML file using safe loader."""
    # Surprisingly, the CSafeLoader does not seem to be used by default.
    # Check whether the CSafeLoader is available and provide a log message
    # if it is not available.
    global _HAS_YAML_LIBRARY

    if _HAS_YAML_LIBRARY is None:
        _HAS_YAML_LIBRARY = hasattr(yaml, 'CSafeLoader')
        if not _HAS_YAML_LIBRARY:
            logger.warning('libyaml was not found! Please install libyaml to'
                           ' speed up loading the model files.')

    if _HAS_YAML_LIBRARY:
        loader = yaml.CSafeLoader(stream)
    else:
        loader = yaml.SafeLoader(stream)
    loader.add_constructor('tag:yaml.org,2002:float', float_constructor)
    return loader.get_data()


class _OrderedEntrySet(object):
    """Ordered set of entity entries.

    This is somewhere between a Set and a Mapping but deliberately does not
    implement any of those interfaces. Iteration does not provide
    (key, value)-items (since the keys are already in the entry values)
    but instead provides just values. The usual set operations such as union
    and intersection are not provided.
    """
    def __init__(self, mapping={}):
        self._dict = OrderedDict(mapping)

    def __contains__(self, key):
        return key in self._dict

    def get(self, key, default=None):
        try:
            return self.__getitem__(key)
        except KeyError:
            return default

    def __getitem__(self, key):
        return self._dict[key]

    def __iter__(self):
        return itervalues(self._dict)

    def __len__(self):
        return len(self._dict)

    def add_entry(self, entry):
        self._dict[entry.id] = entry

    def discard(self, key):
        try:
            del self._dict[key]
        except KeyError:
            pass

    def clear(self):
        self._dict.clear()

    def update(self, it):
        for entry in it:
            self.add_entry(entry)

    def __repr__(self):
        return str('<_OrderedEntrySet {{{}}}>').format(
            ', '.join(repr(x) for x in iter(self)))


class ModelReader(object):
    """Reader of native YAML-based model format.

    The reader can be created from a model YAML file or directly from a
    dict, string or File-like object. Use :meth:`reader_from_path` to read the
    model from a YAML file or directory and use the constructor to read from
    other sources. Any externally referenced file (with ``include``) will be
    read on demand by the parse methods. To read the model fully into memory,
    use the :meth:`create_model` to create a :class:`NativeModel`.
    """

    # Model files to try to open if a directory was specified
    DEFAULT_MODEL = 'model.yaml', 'model.yml'

    def __init__(self, model_from, context=None):
        """Create a model from the specified content.

        Model can be a string, open file, or dictionary.
        """
        if isinstance(model_from, string_types):
            self._model = yaml_load(model_from)
            self._context = context
        elif isinstance(model_from, dict):
            self._context = context
            self._model = model_from
        elif hasattr(model_from, 'read') and callable(model_from.read):
            self._context = context
            self._model = yaml_load(model_from)
        else:
            raise ValueError("Model is of an invalid types")

    @classmethod
    def reader_from_path(cls, path):
        """Create a model from specified path.

        Path can be a directory containing a ``model.yaml`` or ``model.yml``
        file or it can be a path naming the central model file directly.
        """
        context = FilePathContext(path)
        try:
            with context.open('r') as f:
                return ModelReader(f, context)
        except IOError:
            # Try to open the default file
            for filename in cls.DEFAULT_MODEL:
                try:
                    context = FilePathContext(
                        os.path.join(path, filename))
                    with context.open('r') as f:
                        return ModelReader(f, context)
                except:
                    logger.debug('Failed to load model file',
                                 exc_info=True)

        # No model could be loaded
        raise ParseError('No model file could be found ({})'.format(
            ', '.join(cls.DEFAULT_MODEL)))

    @property
    def name(self):
        """Name specified by the model."""
        return self._model.get('name', None)

    @property
    def biomass_reaction(self):
        """Biomass reaction specified by the model."""
        return self._model.get('biomass', None)

    @property
    def extracellular_compartment(self):
        """Extracellular compartment specified by the model.

        Defaults to 'e'.
        """
        return self._model.get('extracellular', 'e')

    @property
    def default_compartment(self):
        """Default compartment specified by the model.

        The compartment that is implied when not specified. In some contexts
        (e.g. for exchange compounds) the extracellular compartment may be
        implied instead. Defaults to 'c'.
        """
        return self._model.get('default_compartment', 'c')

    @property
    def default_flux_limit(self):
        """Default flux limit specified by the model.

        When flux limits on reactions are not specified, this value will be
        used. Flux limit of [0;x] will be implied for irreversible reactions
        and [-x;x] for reversible reactions, where x is this value.
        Defaults to 1000."""
        return self._model.get('default_flux_limit', 1000)

    def parse_compartments(self):
        """Parse compartment information from model.

        Return tuple of: 1) iterator of
        :class:`psamm.datasource.entry.CompartmentEntry`; 2) Set of pairs
        defining the compartment boundaries of the model.
        """

        compartments = OrderedDict()
        boundaries = set()

        if 'compartments' in self._model:
            boundary_map = {}
            for compartment_def in self._model['compartments']:
                compartment_id = convert_to_unicode(compartment_def.get('id'))
                _check_id(compartment_id, 'Compartment')
                if compartment_id in compartments:
                    raise ParseError('Duplicate compartment ID: {}'.format(
                        compartment_id))
                props = dict(compartment_def)
                adjacent_to = props.pop('adjacent_to', None)
                if adjacent_to is not None:
                    if not isinstance(adjacent_to, list):
                        adjacent_to = [adjacent_to]
                    for other in adjacent_to:
                        boundary_map.setdefault(other, set()).add(
                            compartment_id)

                mark = FileMark(self._context, None, None)
                compartment = CompartmentEntry(props, mark)
                compartments[compartment_id] = compartment

            # Check boundaries from boundary_map
            for source, dest_set in iteritems(boundary_map):
                if source not in compartments:
                    raise ParseError(
                        'Invalid compartment {} referenced'
                        ' by compartment {}'.format(
                            source, ', '.join(dest_set)))
                for dest in dest_set:
                    boundaries.add(tuple(sorted((source, dest))))

        return itervalues(compartments), frozenset(boundaries)

    def parse_reactions(self):
        """Yield tuples of reaction ID and reactions defined in the model"""

        # Parse reactions defined in the main model file
        if 'reactions' in self._model:
            for reaction in parse_reaction_list(
                    self._context, self._model['reactions'],
                    self.default_compartment):
                yield reaction

    def has_model_definition(self):
        """Return True when the list of model reactions is set in the model."""
        return 'model' in self._model

    def parse_model(self):
        """Yield reaction IDs of model reactions"""

        if self.has_model_definition():
            for reaction_id in parse_model_group_list(
                    self._context, self._model['model']):
                yield reaction_id

    def parse_limits(self):
        """Yield tuples of reaction ID, lower, and upper bound flux limits"""

        if 'limits' in self._model:
            if not isinstance(self._model['limits'], list):
                raise ParseError('Expected limits to be a list')

            for limit in parse_limits_list(
                    self._context, self._model['limits']):
                yield limit

    def parse_exchange(self):
        """Yield tuples of exchange compounds.

        Each exchange compound is a tuple of compound, reaction ID, lower and
        upper flux limits.
        """

        if 'media' in self._model:
            if 'exchange' in self._model:
                raise ParseError('Both "media" and "exchange" are specified')
            logger.warning(
                'The "media" key is deprecated! Please use "exchange" instead:'
                ' https://psamm.readthedocs.io/en/stable/file_format.html')
            exchange_list = self._model['media']
        else:
            exchange_list = self._model.get('exchange')

        extracellular = self.extracellular_compartment
        if exchange_list is not None:
            if not isinstance(exchange_list, list):
                raise ParseError('Expected "exchange" to be a list')

            for exchange_compound in parse_exchange_list(
                    self._context, exchange_list, extracellular):
                compound, reaction_id, lower, upper = exchange_compound
                if compound.compartment is None:
                    compound = compound.in_compartment(extracellular)
                yield compound, reaction_id, lower, upper

    parse_medium = parse_exchange
    """Yield tuples of exchange compounds.

    .. deprecated:: 0.28
       Use :meth:`parse_exchange` instead.
    """

    def parse_compounds(self):
        """Yield CompoundEntries for defined compounds"""

        if 'compounds' in self._model:
            for compound in parse_compound_list(
                    self._context, self._model['compounds']):
                yield compound

    @property
    def context(self):
        return self._context

    def create_model(self):
        """Return :class:`NativeModel` fully loaded into memory."""
        properties = {
            'name': self.name,
            'biomass': self.biomass_reaction,
            'extracellular': self.extracellular_compartment,
            'default_compartment': self.default_compartment,
            'default_flux_limit': self.default_flux_limit
        }

        if self.context is not None:
            git_version = util.git_try_describe(self.context.basepath)
            properties['version_string'] = git_version

        model = NativeModel(properties)

        # Load compartments into model
        compartment_iter, boundaries = self.parse_compartments()
        for compartment in compartment_iter:
            model.compartments.add_entry(compartment)
        model.compartment_boundaries.update(boundaries)

        # Load compounds into model
        for compound in self.parse_compounds():
            if compound.id in model.compounds:
                existing_entry = model.compounds[compound.id]
                common_props = set(compound.properties).intersection(
                    existing_entry.properties).difference({'id'})
                if len(common_props) > 0:
                    logger.warning(
                        'Compound entry {} at {} overrides already defined'
                        ' properties: {}'.format(
                            compound.id, compound.filemark, common_props))

                properties = dict(compound.properties)
                properties.update(existing_entry.properties)
                compound = CompoundEntry(
                    properties, filemark=compound.filemark)
            model.compounds.add_entry(compound)

        # Load reactions into model
        for reaction in self.parse_reactions():
            if reaction.id in model.reactions:
                existing_entry = model.reactions[reaction.id]
                common_props = set(reaction.properties).intersection(
                    existing_entry.properties).difference({'id'})
                if len(common_props) > 0:
                    logger.warning(
                        'Reaction entry {} at {} overrides already defined'
                        ' properties: {}'.format(
                            reaction.id, reaction.filemark, common_props))

                properties = dict(reaction.properties)
                properties.update(existing_entry.properties)
                reaction = ReactionEntry(
                    properties, filemark=reaction.filemark)
            model.reactions.add_entry(reaction)

        for exchange_def in self.parse_exchange():
            model.exchange[exchange_def[0]] = exchange_def

        for limit in self.parse_limits():
            model.limits[limit[0]] = limit
            dir = model.reactions.get(limit[0]).properties[
                'equation'].direction
            if str(dir) == 'Direction.Forward' or \
                    str(dir) == 'Direction.Right':
                v_min = 0
                v_max = model.default_flux_limit
            elif str(dir) == 'Direction.Both' or str(dir) == 'Direction.Bidir':
                v_min = -model.default_flux_limit
                v_max = model.default_flux_limit
            else:
                logger.warning('Reaction {} has invalid direction: {}'.format(
                    limit[0], dir))
            if (limit[1] is not None and limit[1] < v_min) or \
                    (limit[2] is not None and limit[2] > v_max):
                specify_vmin = limit[1] if limit[1] is not None else \
                    -model.default_flux_limit
                specify_vmax = limit[2] if limit[2] is not None else \
                    -model.default_flux_limit
                logger.warning(
                    "Flux limits for reaction {} override default flux "
                    "limits set by reaction direction (current flux limits: "
                    "[{}, {}]; default limits: [{}, {}]).".format(
                        limit[0], specify_vmin, specify_vmax, v_min, v_max))

        if self.has_model_definition():
            for model_reaction in self.parse_model():
                model.model[model_reaction] = None
        else:
            for reaction in model.reactions:
                model.model[reaction.id] = None

        return model


class NativeModel(object):
    """Represents model in the native format."""

    def __init__(self, properties={}):
        self._properties = dict(properties)
        self._compartments = _OrderedEntrySet()
        self._compartment_boundaries = set()
        self._compounds = _OrderedEntrySet()
        self._reactions = _OrderedEntrySet()
        self._exchange = OrderedDict()
        self._limits = OrderedDict()
        self._model = OrderedDict()

    @property
    def name(self):
        """Return model name property."""
        return self._properties.get('name')

    @name.setter
    def name(self, value):
        self._properties['name'] = value

    @property
    def version_string(self):
        """Return model version string."""
        return self._properties.get('version_string')

    @version_string.setter
    def version_string(self, value):
        self._properties['version_string'] = value

    @property
    def biomass_reaction(self):
        """Return biomass reaction property."""
        return self._properties.get('biomass')

    @biomass_reaction.setter
    def biomass_reaction(self, value):
        self._properties['biomass'] = value

    @property
    def extracellular_compartment(self):
        """Return extracellular compartment property."""
        return self._properties.get('extracellular')

    @extracellular_compartment.setter
    def extracellular_compartment(self, value):
        self._properties['extracellular'] = value

    @property
    def default_compartment(self):
        """Return default compartment property."""
        return self._properties.get('default_compartment')

    @default_compartment.setter
    def default_compartment(self, value):
        self._properties['default_compartment'] = value

    @property
    def default_flux_limit(self):
        """Return default flux limit property."""
        return self._properties.get('default_flux_limit')

    @default_flux_limit.setter
    def default_flux_limit(self, value):
        self._properties['default_flux_limit'] = value

    @property
    def compartments(self):
        """Return compartments entry set."""
        return self._compartments

    @property
    def compartment_boundaries(self):
        """Return set of compartment boundaries."""
        return self._compartment_boundaries

    @property
    def reactions(self):
        """Return reaction entry set."""
        return self._reactions

    @property
    def compounds(self):
        """Return compound entry set."""
        return self._compounds

    @property
    def exchange(self):
        """Return dict of exchange compounds and properties."""
        return self._exchange

    @property
    def limits(self):
        """Return dict of reaction limits."""
        return self._limits

    @property
    def model(self):
        """Return dict of model reactions."""
        return self._model

    def create_metabolic_model(self):
        """Create a :class:`psamm.metabolicmodel.MetabolicModel`."""

        def _translate_compartments(reaction, compartment):
            """Translate compound with missing compartments.

            These compounds will have the specified compartment in the output.
            """
            left = (((c.in_compartment(compartment), v)
                     if c.compartment is None else (c, v))
                    for c, v in reaction.left)
            right = (((c.in_compartment(compartment), v)
                      if c.compartment is None else (c, v))
                     for c, v in reaction.right)
            return Reaction(reaction.direction, left, right)

        # Create metabolic model
        database = DictDatabase()
        for reaction in self.reactions:
            if reaction.equation is not None:
                equation = _translate_compartments(
                    reaction.equation, self.default_compartment)
                database.set_reaction(reaction.id, equation)

        undefined_compartments = set()
        undefined_compounds = set()
        extracellular_compounds = set()
        extracellular = self.extracellular_compartment
        for reaction in database.reactions:
            for compound, _ in database.get_reaction_values(reaction):
                if compound.name not in self.compounds:
                    undefined_compounds.add(compound.name)
                if compound.compartment == extracellular:
                    extracellular_compounds.add(compound.name)
                if compound.compartment not in self.compartments:
                    undefined_compartments.add(compound.compartment)

        for compartment in sorted(undefined_compartments):
            logger.warning(
                'The compartment {} was not defined in the list'
                ' of compartments'.format(compartment))

        for compound in sorted(undefined_compounds):
            logger.warning(
                'The compound {} was not defined in the list'
                ' of compounds'.format(compound))

        exchange_compounds = set()
        for exchange_compound in self.exchange:
            if exchange_compound.compartment == extracellular:
                exchange_compounds.add(exchange_compound.name)

        for compound in sorted(extracellular_compounds - exchange_compounds):
            logger.warning(
                'The compound {} was in the extracellular compartment'
                ' but not defined in the exchange compounds'.format(compound))
        for compound in sorted(exchange_compounds - extracellular_compounds):
            logger.warning(
                'The compound {} was defined in the exchange compounds but'
                ' is not in the extracellular compartment'.format(compound))

        model_definition = None
        if len(self.model) > 0:
            model_definition = self.model

        for reaction in self.reactions:
            if reaction.equation is None:
                logger.warning(
                    'Reaction {} has no reaction equation'.format(reaction.id))
                del model_definition[reaction.id]

        return MetabolicModel.load_model(
            database, model_definition, itervalues(self.exchange),
            itervalues(self.limits), v_max=self.default_flux_limit)

    def __repr__(self):
        return str('<{} name={!r}>'.format(self.__class__.__name__, self.name))


def _check_id(entity, entity_type):
    """Check whether the ID is valid.

    First check if the ID is missing, and then check if it is a qualified
    string type, finally check if the string is empty. For all checks, it
    would raise a ParseError with the corresponding message.

    Args:
        entity: a string type object to be checked.
        entity_type: a string that shows the type of entities to check, usually
            `Compound` or 'Reaction'.
    """

    if entity is None:
        raise ParseError('{} ID missing'.format(entity_type))
    elif not isinstance(entity, string_types):
        msg = '{} ID must be a string, id was {}.'.format(entity_type, entity)
        if isinstance(entity, bool):
            msg += (' You may have accidentally used an ID value that YAML'
                    ' interprets as a boolean, such as "yes", "no", "on",'
                    ' "off", "true" or "false". To use this ID, you have to'
                    ' quote it with single or double quotes')
        raise ParseError(msg)
    elif len(entity) == 0:
        raise ParseError('{} ID must not be empty'.format(entity_type))


def parse_compound(compound_def, context=None):
    """Parse a structured compound definition as obtained from a YAML file

    Returns a CompoundEntry."""
    compound_id = convert_to_unicode(compound_def.get('id'))
    _check_id(compound_id, 'Compound')

    compound_props = {}
    for key, value in compound_def.items():
        if isinstance(value, str):
            compound_props[key] = convert_to_unicode(value)
        else:
            compound_props[key] = value

    mark = FileMark(context, None, None)
    return CompoundEntry(compound_props, mark)


def parse_compound_list(path, compounds):
    """Parse a structured list of compounds as obtained from a YAML file

    Yields CompoundEntries. Path can be given as a string or a context.
    """
    context = FilePathContext(path)

    for compound_def in compounds:
        compound_dict = dict(compound_def)
        if 'include' in compound_dict:
            file_format = compound_dict.get('format')
            include_context = context.resolve(compound_dict['include'])
            for compound in parse_compound_file(include_context, file_format):
                yield compound
        else:
            yield parse_compound(compound_dict, context)


def parse_compound_table_file(path, f):
    """Parse a tab-separated file containing compound IDs and properties

    The compound properties are parsed according to the header which specifies
    which property is contained in each column.
    """

    context = FilePathContext(path)

    for i, row in enumerate(csv.DictReader(f, delimiter=str('\t'))):
        if text_type('id') not in row or \
                text_type(convert_to_unicode(row['id'].strip())) == '':
            raise ParseError('Expected `id` column in table')

        props = {key: text_type(convert_to_unicode(value))
                 for key, value in iteritems(row)
                 if text_type(convert_to_unicode(value)) != ''}
        if 'charge' in props:
            props['charge'] = int(props['charge'])

        mark = FileMark(context, i + 2, None)
        yield CompoundEntry(props, mark)


def parse_compound_yaml_file(path, f):
    """Parse a file as a YAML-format list of compounds

    Path can be given as a string or a context.
    """
    return parse_compound_list(path, yaml_load(f))


def resolve_format(format, path):
    """Looks at a file's extension and format (if any) and returns format.
    """
    if format is None:
        if (re.match(r'.+\.(yml|yaml)$', path)):
            return 'yaml'
        elif (re.match(r'.+\.tsv$', path)):
            return 'tsv'
    else:
        return format.lower()


def parse_compound_file(path, format):
    """Open and parse reaction file based on file extension or given format

    Path can be given as a string or a context.
    """

    context = FilePathContext(path)

    # YAML files do not need to explicitly specify format
    format = resolve_format(format, context.filepath)
    if format == 'yaml':
        logger.debug('Parsing compound file {} as YAML'.format(
            context.filepath))
        with context.open('r') as f:
            for compound in parse_compound_yaml_file(context, f):
                yield compound
    elif format == 'modelseed':
        logger.debug('Parsing compound file {} as ModelSEED TSV'.format(
            context.filepath))
        with context.open('r') as f:
            for compound in modelseed.parse_compound_file(f, context):
                yield compound
    elif format == 'tsv':
        logger.debug('Parsing compound file {} as TSV'.format(
            context.filepath))
        with context.open('r') as f:
            for compound in parse_compound_table_file(context, f):
                yield compound
    else:
        raise ParseError('Unable to detect format of compound file {}'.format(
            context.filepath))


def parse_reaction_equation_string(equation, default_compartment):
    """Parse a string representation of a reaction equation.

    Converts undefined compartments to the default compartment.
    """
    def _translate_compartments(reaction, compartment):
        """Translate compound with missing compartments.

        These compounds will have the specified compartment in the output.
        """
        left = (((c.in_compartment(compartment), v)
                 if c.compartment is None else (c, v))
                for c, v in reaction.left)
        right = (((c.in_compartment(compartment), v)
                  if c.compartment is None else (c, v))
                 for c, v in reaction.right)
        return Reaction(reaction.direction, left, right)
    eq = _REACTION_PARSER.parse(equation).normalized()
    return _translate_compartments(eq, default_compartment)


def parse_reaction_equation(equation_def, default_compartment):
    """Parse a structured reaction equation as obtained from a YAML file

    Returns a Reaction.
    """
    def parse_compound_list(compound_list, compartment):
        """Parse a list of reactants or metabolites"""
        for compound_def in compound_list:
            compound_id = convert_to_unicode(compound_def.get('id'))
            _check_id(compound_id, 'Compound')

            value = compound_def.get('value')
            if value is None:
                raise ParseError('Missing value for compound {}'.format(
                    compound_id))

            compound_compartment = compound_def.get('compartment')
            if compound_compartment is None:
                compound_compartment = compartment

            compound = Compound(compound_id, compartment=compound_compartment)
            yield compound, value

    if isinstance(equation_def, string_types):
        return parse_reaction_equation_string(
            equation_def, default_compartment)
    else:
        compartment = equation_def.get('compartment', default_compartment)
        reversible = bool(equation_def.get('reversible', True))
        left = equation_def.get('left', [])
        right = equation_def.get('right', [])
        if len(left) == 0 and len(right) == 0:
            raise ParseError('Reaction values are missing')

        return Reaction(Direction.Both if reversible else Direction.Forward,
                        parse_compound_list(left, compartment),
                        parse_compound_list(right, compartment))


def parse_reaction(reaction_def, default_compartment, context=None):
    """Parse a structured reaction definition as obtained from a YAML file

    Returns a ReactionEntry.
    """
    reaction_id = convert_to_unicode(reaction_def.get('id'))
    _check_id(reaction_id, 'Reaction')

    reaction_props = {}
    for key, value in reaction_def.items():
        if isinstance(value, str):
            reaction_props[key] = convert_to_unicode(value)
        else:
            reaction_props[key] = value

    if 'genes' in reaction_def.keys():
        if isinstance(reaction_def['genes'], list):
            reaction_props['genes'] = [convert_to_unicode(gene)
                                       for gene in reaction_def['genes']]

    # Parse reaction equation
    if 'equation' in reaction_def.keys():
        reaction_props['equation'] = parse_reaction_equation(
            reaction_def['equation'], default_compartment)

    mark = FileMark(context, None, None)
    return ReactionEntry(reaction_props, mark)


def parse_reaction_list(path, reactions, default_compartment=None):
    """Parse a structured list of reactions as obtained from a YAML file

    Yields tuples of reaction ID and reaction object. Path can be given as a
    string or a context.
    """
    context = FilePathContext(path)

    for reaction_def in reactions:
        reaction_dict = dict(reaction_def)
        if 'include' in reaction_dict:
            include_context = context.resolve(reaction_dict['include'])
            for reaction in parse_reaction_file(
                    include_context, default_compartment):
                yield reaction
        else:
            yield parse_reaction(reaction_dict, default_compartment, context)


def parse_reaction_yaml_file(path, f, default_compartment):
    """Parse a file as a YAML-format list of reactions

    Path can be given as a string or a context.
    """

    return parse_reaction_list(path, yaml_load(f), default_compartment)


def parse_reaction_table_file(path, f, default_compartment):
    """Parse a tab-separated file containing reaction IDs and properties

    The reaction properties are parsed according to the header which specifies
    which property is contained in each column.
    """

    context = FilePathContext(path)

    for lineno, row in enumerate(csv.DictReader(f, delimiter=str('\t'))):
        if text_type('id') not in row or \
                text_type(convert_to_unicode(row['id'].strip())) == '':
            raise ParseError('Expected `id` column in table')

        props = {key: convert_to_unicode(value)
                 for key, value in iteritems(row)
                 if convert_to_unicode(value) != ''}

        if 'equation' in props:
            props['equation'] = parse_reaction_equation_string(
                props['equation'], default_compartment)

        mark = FileMark(context, lineno + 2, 0)
        yield ReactionEntry(props, mark)


def parse_reaction_file(path, default_compartment=None):
    """Open and parse reaction file based on file extension

    Path can be given as a string or a context.
    """

    context = FilePathContext(path)

    format = resolve_format(None, context.filepath)
    if format == 'tsv':
        logger.debug('Parsing reaction file {} as TSV'.format(
            context.filepath))
        with context.open('r') as f:
            for reaction in parse_reaction_table_file(
                    context, f, default_compartment):
                yield reaction
    elif format == 'yaml':
        logger.debug('Parsing reaction file {} as YAML'.format(
            context.filepath))
        with context.open('r') as f:
            for reaction in parse_reaction_yaml_file(
                    context, f, default_compartment):
                yield reaction
    else:
        raise ParseError('Unable to detect format of reaction file {}'.format(
            context.filepath))


def get_limits(compound_def):
    if ('fixed' in compound_def and
            ('lower' not in compound_def and 'upper'not in compound_def)):
        fixed = compound_def['fixed']
        lower = fixed
        upper = fixed
    elif ('fixed' in compound_def and
            ('lower' in compound_def or 'upper' in compound_def)):
        raise ParseError('Cannot use fixed and a lower or upper bound')
    else:
        lower = compound_def.get('lower', None)
        upper = compound_def.get('upper', None)
    return lower, upper


def parse_exchange(exchange_def, default_compartment):
    """Parse a structured exchange definition as obtained from a YAML file.

    Returns in iterator of compound, reaction, lower and upper bounds.
    """

    default_compartment = exchange_def.get('compartment', default_compartment)

    for compound_def in exchange_def.get('compounds', []):
        compartment = compound_def.get('compartment', default_compartment)
        compound = Compound(convert_to_unicode(compound_def['id']),
                            compartment=compartment)
        reaction = compound_def.get('reaction')
        if reaction:
            reaction = convert_to_unicode(reaction)
        lower, upper = get_limits(compound_def)
        yield compound, reaction, lower, upper


parse_medium = parse_exchange
"""Parse a structured exchange definition as obtained from a YAML file.

.. deprecated:: 0.28
   Use :func:`parse_exchange` instead.
"""


def parse_exchange_list(path, exchange, default_compartment):
    """Parse a structured exchange list as obtained from a YAML file.

    Yields tuples of compound, reaction ID, lower and upper flux bounds. Path
    can be given as a string or a context.
    """

    context = FilePathContext(path)

    for exchange_def in exchange:
        if 'include' in exchange_def:
            include_context = context.resolve(exchange_def['include'])
            for exchange_compound in parse_exchange_file(
                    include_context, default_compartment):
                yield exchange_compound
        else:
            for exchange_compound in parse_exchange(
                    exchange_def, default_compartment):
                yield exchange_compound


parse_medium_list = parse_exchange_list
"""Parse a structured exchange list as obtained from a YAML file.

.. deprecated:: 0.28
   Use :func:`parse_exchange_list` instead.
"""


def parse_exchange_yaml_file(path, f, default_compartment):
    """Parse a file as a YAML-format exchange definition.

    Path can be given as a string or a context.
    """

    return parse_exchange(yaml_load(f), default_compartment)


parse_medium_yaml_file = parse_exchange_yaml_file
"""Parse a file as a YAML-format exchange definition.

.. deprecated:: 0.28
   Use :func:`parse_exchange_yaml_file` instead.
"""


def parse_exchange_table_file(f):
    """Parse a space-separated file containing exchange compound flux limits.

    The first two columns contain compound IDs and compartment while the
    third column contains the lower flux limits. The fourth column is
    optional and contains the upper flux limit.
    """

    for line in f:
        line, _, comment = convert_to_unicode(line).partition('#')
        line = line.strip()
        if line == '':
            continue

        # A line can specify lower limit only (useful for
        # medium compounds), or both lower and upper limit.
        fields = line.split(None)
        if len(fields) < 2 or len(fields) > 4:
            raise ParseError('Malformed compound limit: {}'.format(fields))

        # Extend to four values and unpack
        fields.extend(['-']*(4-len(fields)))
        compound_id, compartment, lower, upper = fields

        compound = Compound(convert_to_unicode(compound_id), compartment)
        lower = float(lower) if lower != '-' else None
        upper = float(upper) if upper != '-' else None

        yield compound, None, lower, upper


parse_medium_table_file = parse_exchange_table_file
"""Parse a space-separated file containing exchange compound flux limits.

.. deprecated:: 0.28
   Use :func:`parse_exchange_table_file` instead.
"""


def parse_exchange_file(path, default_compartment):
    """Parse a file as a list of exchange compounds with flux limits.

    The file format is detected and the file is parsed accordingly. Path can
    be given as a string or a context.
    """

    context = FilePathContext(path)

    format = resolve_format(None, context.filepath)
    if format == 'tsv':
        logger.debug('Parsing exchange file {} as TSV'.format(
            context.filepath))
        with context.open('r') as f:
            for entry in parse_exchange_table_file(f):
                yield entry
    elif format == 'yaml':
        logger.debug('Parsing exchange file {} as YAML'.format(
            context.filepath))
        with context.open('r') as f:
            for entry in parse_exchange_yaml_file(
                    context, f, default_compartment):
                yield entry
    else:
        raise ParseError('Unable to detect format of exchange file {}'.format(
            context.filepath))


parse_medium_file = parse_exchange_file
"""Parse a file as a list of exchange compounds with flux limits.

.. deprecated:: 0.28
   Use :func:`parse_exchange_file` instead.
"""


def parse_limit(limit_def):
    """Parse a structured flux limit definition as obtained from a YAML file

    Returns a tuple of reaction, lower and upper bound.
    """

    lower, upper = get_limits(limit_def)
    reaction = convert_to_unicode(limit_def.get('reaction'))

    return reaction, lower, upper


def parse_limits_list(path, limits):
    """Parse a structured list of flux limits as obtained from a YAML file

    Yields tuples of reaction ID, lower and upper flux bounds. Path can be
    given as a string or a context.
    """

    context = FilePathContext(path)

    for limit_def in limits:
        if 'include' in limit_def:
            include_context = context.resolve(limit_def['include'])
            for limit in parse_limits_file(include_context):
                yield limit
        else:
            yield parse_limit(limit_def)


def parse_limits_table_file(f):
    """Parse a space-separated file containing reaction flux limits

    The first column contains reaction IDs while the second column contains
    the lower flux limits. The third column is optional and contains the
    upper flux limit.
    """

    for line in f:
        line, _, comment = convert_to_unicode(line).partition('#')
        line = line.strip()
        if line == '':
            continue

        # A line can specify lower limit only (useful for
        # exchange reactions), or both lower and upper limit.
        fields = line.split(None)
        if len(fields) < 1 or len(fields) > 3:
            raise ParseError('Malformed reaction limit: {}'.format(fields))

        # Extend to three values and unpack
        fields.extend(['-']*(3-len(fields)))
        reaction_id, lower, upper = fields

        lower = float(lower) if lower != '-' else None
        upper = float(upper) if upper != '-' else None

        yield convert_to_unicode(reaction_id), lower, upper


def parse_limits_yaml_file(path, f):
    """Parse a file as a YAML-format flux limits definition

    Path can be given as a string or a context.
    """

    return parse_limits_list(path, yaml_load(f))


def parse_limits_file(path):
    """Parse a file as a list of reaction flux limits

    The file format is detected and the file is parsed accordingly. Path can
    be given as a string or a context.
    """

    context = FilePathContext(path)

    format = resolve_format(None, context.filepath)
    if format == 'tsv':
        logger.debug('Parsing limits file {} as TSV'.format(
            context.filepath))
        with context.open('r') as f:
            for limit in parse_limits_table_file(f):
                yield limit
    elif format == 'yaml':
        logger.debug('Parsing limits file {} as YAML'.format(
            context.filepath))
        with context.open('r') as f:
            for limit in parse_limits_yaml_file(context, f):
                yield limit
    else:
        raise ParseError('Unable to detect format of limits file {}'.format(
            context.filepath))


def parse_model_group(path, group):
    """Parse a structured model group as obtained from a YAML file

    Path can be given as a string or a context.
    """

    context = FilePathContext(path)

    for reaction_id in group.get('reactions', []):
        yield convert_to_unicode(reaction_id)

    # Parse subgroups
    for reaction_id in parse_model_group_list(
            context, group.get('groups', [])):
        yield reaction_id


def parse_model_group_list(path, groups):
    """Parse a structured list of model groups as obtained from a YAML file

    Yields reaction IDs. Path can be given as a string or a context.
    """

    context = FilePathContext(path)
    for model_group in groups:
        if 'include' in model_group:
            include_context = context.resolve(model_group['include'])
            for reaction_id in parse_model_file(include_context):
                yield reaction_id
        else:
            for reaction_id in parse_model_group(context, model_group):
                yield reaction_id


def parse_model_yaml_file(path, f):
    """Parse a file as a YAML-format list of model reaction groups

    Path can be given as a string or a context.
    """
    return parse_model_group_list(path, yaml_load(f))


def parse_model_table_file(path, f):
    """Parse a file as a list of model reactions

    Yields reactions IDs. Path can be given as a string or a context.
    """

    for line in f:
        line, _, comment = convert_to_unicode(line).partition('#')
        line = line.strip()
        if line == '':
            continue

        yield line


def parse_model_file(path):
    """Parse a file as a list of model reactions

    The file format is detected and the file is parsed accordinly. The file is
    specified as a file path that will be opened for reading. Path can be given
    as a string or a context.
    """

    context = FilePathContext(path)

    format = resolve_format(None, context.filepath)
    if format == 'tsv':
        logger.debug('Parsing model file {} as TSV'.format(context.filepath))
        with context.open('r') as f:
            for reaction_id in parse_model_table_file(context, f):
                yield reaction_id
    elif format == 'yaml':
        logger.debug('Parsing model file {} as YAML'.format(context.filepath))
        with context.open('r') as f:
            for reaction_id in parse_model_yaml_file(context, f):
                yield reaction_id


# Threshold for converting reactions into dictionary representation.
_MAX_REACTION_LENGTH = 10


# Create wrapper representer for text_type for Py2/3 compatibility.
if PY3:
    def _represent_text_type(dumper, data):
        return dumper.represent_str(data)
else:
    def _represent_text_type(dumper, data):
        return dumper.represent_unicode(data)


# Define custom dict representers for YAML
# This allows reading/writing Python OrderedDicts in the correct order.
# See: https://stackoverflow.com/questions/5121931/in-python-how-can-you-load-yaml-mappings-as-ordereddicts  # noqa
def _dict_representer(dumper, data):
    return dumper.represent_dict(iteritems(data))


def _set_representer(dumper, data):
    return dumper.represent_list(iter(data))


def _boolean_expression_representer(dumper, data):
    return _represent_text_type(dumper, text_type(data))


def _reaction_representer(dumper, data):
    """Generate a parsable reaction representation to the YAML parser.

    Check the number of compounds in the reaction, if it is larger than 10,
    then transform the reaction data into a list of directories with all
    attributes in the reaction; otherwise, just return the text_type format
    of the reaction data.
    """
    if len(data.compounds) > _MAX_REACTION_LENGTH:
        def dict_make(compounds):
            for compound, value in compounds:
                yield OrderedDict([
                    ('id', text_type(compound.name)),
                    ('compartment', compound.compartment),
                    ('value', value)])

        left = list(dict_make(data.left))
        right = list(dict_make(data.right))

        direction = data.direction == Direction.Both

        reaction = OrderedDict()
        reaction['reversible'] = direction
        if data.direction == Direction.Reverse:
            reaction['left'] = right
            reaction['right'] = left
        else:
            reaction['left'] = left
            reaction['right'] = right

        return dumper.represent_data(reaction)
    else:
        return _represent_text_type(dumper, text_type(data))


def _formula_representer(dumper, data):
    return _represent_text_type(dumper, text_type(data))


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


class ModelWriter(object):
    """Writer for native (YAML) format."""

    def __init__(self):
        self._yaml_args = {
            'default_flow_style': False,
            'encoding': 'utf-8',
            'allow_unicode': True,
            'width': 79
        }

    def _dump(self, stream, data):
        if hasattr(yaml, 'CSafeDumper'):
            dumper = yaml.CSafeDumper(stream, **self._yaml_args)
        else:
            dumper = yaml.SafeDumper(stream, **self._yaml_args)

        dumper.add_representer(OrderedDict, _dict_representer)
        dumper.add_representer(set, _set_representer)
        dumper.add_representer(frozenset, _set_representer)
        dumper.add_representer(
            boolean.Expression, _boolean_expression_representer)
        dumper.add_representer(Reaction, _reaction_representer)
        dumper.add_representer(Formula, _formula_representer)
        dumper.add_representer(Decimal, _decimal_representer)

        dumper.ignore_aliases = lambda *args: True

        try:
            dumper.open()
            dumper.represent(data)
            dumper.close()
        finally:
            dumper.dispose()

    def convert_compartment_entry(self, compartment, adjacencies):
        """Convert compartment entry to YAML dict.

        Args:
            compartment: :class:`psamm.datasource.entry.CompartmentEntry`.
            adjacencies: Sequence of IDs or a single ID of adjacent
                compartments (or None).
        """
        d = OrderedDict()
        d['id'] = compartment.id
        if adjacencies is not None:
            d['adjacent_to'] = adjacencies

        order = {key: i for i, key in enumerate(['name'])}
        prop_keys = set(compartment.properties)
        for prop in sorted(prop_keys,
                           key=lambda x: (order.get(x, 1000), x)):
            if compartment.properties[prop] is not None:
                d[prop] = compartment.properties[prop]

        return d

    def convert_compound_entry(self, compound):
        """Convert compound entry to YAML dict."""
        d = OrderedDict()
        d['id'] = compound.id

        order = {
            key: i for i, key in enumerate(
                ['name', 'formula', 'formula_neutral', 'charge', 'kegg',
                 'cas'])}
        prop_keys = (
            set(compound.properties) - {'boundary', 'compartment'})
        for prop in sorted(prop_keys,
                           key=lambda x: (order.get(x, 1000), x)):
            if compound.properties[prop] is not None:
                d[prop] = compound.properties[prop]

        return d

    def convert_reaction_entry(self, reaction):
        """Convert reaction entry to YAML dict."""
        d = OrderedDict()
        d['id'] = reaction.id

        def is_equation_valid(equation):
            # If the equation is a Reaction object, it must have non-zero
            # number of compounds.
            return (equation is not None and (
                    not isinstance(equation, Reaction) or
                    len(equation.compounds) > 0))

        order = {
            key: i for i, key in enumerate(
                ['name', 'genes', 'equation', 'subsystem', 'ec'])}
        prop_keys = (set(reaction.properties) -
                     {'lower_flux', 'upper_flux', 'reversible'})
        for prop in sorted(prop_keys, key=lambda x: (order.get(x, 1000), x)):
            if reaction.properties[prop] is None:
                continue
            d[prop] = reaction.properties[prop]
            if prop == 'equation' and not is_equation_valid(d[prop]):
                del d[prop]

        return d

    def _write_entries(self, stream, entries, converter, properties=None):
        """Write iterable of entries as YAML object to stream.

        Args:
            stream: File-like object.
            entries: Iterable of entries.
            converter: Conversion function from entry to YAML object.
            properties: Set of compartment properties to output (or None to
                output all).
        """
        def iter_entries():
            for c in entries:
                entry = converter(c)
                if entry is None:
                    continue
                if properties is not None:
                    entry = OrderedDict(
                        (key, value) for key, value in iteritems(entry)
                        if key == 'id' or key in properties)
                yield entry

        self._dump(stream, list(iter_entries()))

    def write_compartments(self, stream, compartments, adjacencies,
                           properties=None):
        """Write iterable of compartments as YAML object to stream.

        Args:
            stream: File-like object.
            compartments: Iterable of compartment entries.
            adjacencies: Dictionary mapping IDs to adjacent compartment IDs.
            properties: Set of compartment properties to output (or None to
                output all).
        """
        def convert(entry):
            return self.convert_compartment_entry(
                entry, adjacencies.get(entry.id))

        self._write_entries(stream, compartments, convert, properties)

    def write_compounds(self, stream, compounds, properties=None):
        """Write iterable of compounds as YAML object to stream.

        Args:
            stream: File-like object.
            compounds: Iterable of compound entries.
            properties: Set of compound properties to output (or None to output
                all).
        """
        self._write_entries(
            stream, compounds, self.convert_compound_entry, properties)

    def write_reactions(self, stream, reactions, properties=None):
        """Write iterable of reactions as YAML object to stream.

        Args:
            stream: File-like object.
            compounds: Iterable of reaction entries.
            properties: Set of reaction properties to output (or None to output
                all).
        """
        self._write_entries(
            stream, reactions, self.convert_reaction_entry, properties)


def reaction_signature(eq, direction=False, stoichiometry=False):
    """Return unique signature object for :class:`Reaction`.

    Signature objects are hashable, and compare equal only if the reactions
    are considered the same according to the specified rules.

    Args:
        direction: Include reaction directionality when considering equality.
        stoichiometry: Include stoichiometry when considering equality.
    """
    def compounds_sig(compounds):
        if stoichiometry:
            return tuple(sorted(compounds))
        else:
            return tuple(sorted(compound for compound, _ in compounds))

    left = compounds_sig(eq.left)
    right = compounds_sig(eq.right)

    if left < right:
        reaction_sig = left, right
        direction_sig = eq.direction
    else:
        reaction_sig = right, left
        direction_sig = eq.direction.flipped()

    if direction:
        return reaction_sig, direction_sig
    return reaction_sig
