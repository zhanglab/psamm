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
# Copyright 2014-2016  Jon Lund Steffensen <jon_steffensen@uri.edu>

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
from six import string_types, text_type, iteritems, PY3
from decimal import Decimal

from ..reaction import Reaction, Compound, Direction
from ..expression import boolean
from ..formula import Formula
from ..metabolicmodel import MetabolicModel
from ..database import DictDatabase
from .context import FilePathContext, FileMark
from .entry import (DictCompoundEntry as CompoundEntry,
                    DictReactionEntry as ReactionEntry)
from .reaction import ReactionParser
from . import modelseed

# Module-level logging
logger = logging.getLogger(__name__)

# Model files to try to open if a directory was specified
DEFAULT_MODEL = ('model.yaml', 'model.yml')

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


def whendefined(func, value):
    """Apply func to value if value is not None"""
    return func(value) if value is not None else None


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


class NativeModel(object):
    """Represents a model specified using the native data formats

    The model is created from a model file or from a directory containing a
    model file using the default file name (model.yaml or model.yml). This file
    can specify the model fully or refer to other files within the same
    directory subtree that specifies part of the model.
    """

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
    def load_model_from_path(cls, path):
        """Create a model from specified path."""
        context = FilePathContext(path)
        try:
            with open(context.filepath, 'r') as f:
                return NativeModel(f, context)
        except IOError:
            # Try to open the default file
            for filename in DEFAULT_MODEL:
                try:
                    context = FilePathContext(
                        os.path.join(path, filename))
                    with open(context.filepath, 'r') as f:
                        return NativeModel(f, context)
                except:
                    logger.debug('Failed to load model file',
                                 exc_info=True)

        # No model could be loaded
        raise ParseError('No model file could be found ({})'.format(
            ', '.join(DEFAULT_MODEL)))

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
        (e.g. medium) the extracellular compartment may be implied instead.
        Defaults to 'c'.
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

    def parse_medium(self):
        """Yield tuples of medium compounds.

        Each medium compound is a tuple of compound, reaction ID, lower and
        upper flux limits.
        """

        extracellular = self.extracellular_compartment
        if 'media' in self._model:
            if not isinstance(self._model['media'], list):
                raise ParseError('Expected media to be a list')

            for medium_compound in parse_medium_list(
                    self._context, self._model['media'], extracellular):
                compound, reaction_id, lower, upper = medium_compound
                if compound.compartment is None:
                    compound = compound.in_compartment(extracellular)
                yield compound, reaction_id, lower, upper

    def parse_compounds(self):
        """Yield CompoundEntries for defined compounds"""

        if 'compounds' in self._model:
            for compound in parse_compound_list(
                    self._context, self._model['compounds']):
                yield compound

    def create_metabolic_model(self):
        """Create a :class:`psamm.metabolicmodel.MetabolicModel`."""

        # Create metabolic model
        database = DictDatabase()
        for reaction in self.parse_reactions():
            if reaction.equation is not None:
                database.set_reaction(reaction.id, reaction.equation)

        # Warn about undefined compounds
        compounds = set()
        for compound in self.parse_compounds():
            compounds.add(compound.id)

        undefined_compounds = set()
        extracellular_compounds = set()
        extracellular = self.extracellular_compartment
        for reaction in database.reactions:
            for compound, _ in database.get_reaction_values(reaction):
                if compound.name not in compounds:
                    undefined_compounds.add(compound.name)
                if compound.compartment == extracellular:
                    extracellular_compounds.add(compound.name)

        for compound in sorted(undefined_compounds):
            logger.warning(
                'The compound {} was not defined in the list'
                ' of compounds'.format(compound))

        medium_compounds = set()
        for medium_compound in self.parse_medium():
            if medium_compound[0].compartment == extracellular:
                medium_compounds.add(medium_compound[0].name)

        for compound in sorted(extracellular_compounds - medium_compounds):
            logger.warning(
                'The compound {} was in the extracellular compartment'
                ' but not defined in the medium'.format(compound))
        for compound in sorted(medium_compounds - extracellular_compounds):
            logger.warning(
                'The compound {} was defined in the medium but'
                ' is not in the extracellular compartment'.format(compound))

        model_definition = None
        if self.has_model_definition():
            model_definition = self.parse_model()

        return MetabolicModel.load_model(
            database, model_definition, self.parse_medium(),
            self.parse_limits(), v_max=self.default_flux_limit)

    @property
    def context(self):
        """Model file loading context (or None, if undefined)."""
        return self._context


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

    compound_id = compound_def.get('id')
    _check_id(compound_id, 'Compound')

    mark = FileMark(context, None, None)
    return CompoundEntry(compound_def, mark)


def parse_compound_list(path, compounds):
    """Parse a structured list of compounds as obtained from a YAML file

    Yields CompoundEntries. Path can be given as a string or a context.
    """

    context = FilePathContext(path)

    for compound_def in compounds:
        if 'include' in compound_def:
            file_format = compound_def.get('format')
            include_context = context.resolve(compound_def['include'])
            for compound in parse_compound_file(include_context, file_format):
                yield compound
        else:
            yield parse_compound(compound_def, context)


def parse_compound_table_file(path, f):
    """Parse a tab-separated file containing compound IDs and properties

    The compound properties are parsed according to the header which specifies
    which property is contained in each column.
    """

    context = FilePathContext(path)

    for i, row in enumerate(csv.DictReader(f, delimiter=str('\t'))):
        if 'id' not in row or row['id'].strip() == '':
            raise ParseError('Expected `id` column in table')

        props = {key: value for key, value in iteritems(row) if value != ''}

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

    def parse_compound_list(l, compartment):
        """Parse a list of reactants or metabolites"""
        for compound_def in l:
            compound_id = compound_def.get('id')
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

    reaction_id = reaction_def.get('id')
    _check_id(reaction_id, 'Reaction')

    reaction_props = dict(reaction_def)

    # Parse reaction equation
    if 'equation' in reaction_def:
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
        if 'include' in reaction_def:
            include_context = context.resolve(reaction_def['include'])
            for reaction in parse_reaction_file(
                    include_context, default_compartment):
                yield reaction
        else:
            yield parse_reaction(reaction_def, default_compartment, context)


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
        if 'id' not in row or row['id'].strip() == '':
            raise ParseError('Expected `id` column in table')

        props = {key: value for key, value in iteritems(row) if value != ''}

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
            ('lower'in compound_def or 'upper' in compound_def)):
        raise ParseError('Cannot use fixed and a lower or upper bound')
    else:
        lower = compound_def.get('lower', None)
        upper = compound_def.get('upper', None)
    return lower, upper


def parse_medium(medium_def, default_compartment):
    """Parse a structured medium definition as obtained from a YAML file

    Returns in iterator of compound, reaction, lower and upper bounds.
    """

    default_compartment = medium_def.get('compartment', default_compartment)

    for compound_def in medium_def.get('compounds', []):
        compartment = compound_def.get('compartment', default_compartment)
        compound = Compound(compound_def['id'], compartment=compartment)
        reaction = compound_def.get('reaction')
        lower, upper = get_limits(compound_def)
        yield compound, reaction, lower, upper


def parse_medium_list(path, medium, default_compartment):
    """Parse a structured medium list as obtained from a YAML file.

    Yields tuples of compound, reaction ID, lower and upper flux bounds. Path
    can be given as a string or a context.
    """

    context = FilePathContext(path)

    for medium_def in medium:
        if 'include' in medium_def:
            include_context = context.resolve(medium_def['include'])
            for medium_compound in parse_medium_file(
                    include_context, default_compartment):
                yield medium_compound
        else:
            for medium_compound in parse_medium(
                    medium_def, default_compartment):
                yield medium_compound


def parse_medium_yaml_file(path, f, default_compartment):
    """Parse a file as a YAML-format medium definition

    Path can be given as a string or a context.
    """

    return parse_medium(yaml_load(f), default_compartment)


def parse_medium_table_file(f):
    """Parse a space-separated file containing medium compound flux limits

    The first two columns contain compound IDs and compartment while the
    third column contains the lower flux limits. The fourth column is
    optional and contains the upper flux limit.
    """

    for line in f:
        line, _, comment = line.partition('#')
        line = line.strip()
        if line == '':
            continue

        # A line can specify lower limit only (useful for
        # exchange reactions), or both lower and upper limit.
        fields = line.split(None)
        if len(fields) < 2 or len(fields) > 4:
            raise ParseError('Malformed compound limit: {}'.format(fields))

        # Extend to four values and unpack
        fields.extend(['-']*(4-len(fields)))
        compound_id, compartment, lower, upper = fields

        compound = Compound(compound_id, compartment)
        lower = float(lower) if lower != '-' else None
        upper = float(upper) if upper != '-' else None

        yield compound, None, lower, upper


def parse_medium_file(path, default_compartment):
    """Parse a file as a list of medium compounds with flux limits

    The file format is detected and the file is parsed accordingly. Path can
    be given as a string or a context.
    """

    context = FilePathContext(path)

    format = resolve_format(None, context.filepath)
    if format == 'tsv':
        logger.debug('Parsing medium file {} as TSV'.format(
            context.filepath))
        with context.open('r') as f:
            for entry in parse_medium_table_file(f):
                yield entry
    elif format == 'yaml':
        logger.debug('Parsing medium file {} as YAML'.format(
            context.filepath))
        with context.open('r') as f:
            for entry in parse_medium_yaml_file(
                    context, f, default_compartment):
                yield entry
    else:
        raise ParseError('Unable to detect format of medium file {}'.format(
            context.filepath))


def parse_limit(limit_def):
    """Parse a structured flux limit definition as obtained from a YAML file

    Returns a tuple of reaction, lower and upper bound.
    """

    lower, upper = get_limits(limit_def)
    reaction = limit_def.get('reaction')

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
        line, _, comment = line.partition('#')
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

        yield reaction_id, lower, upper


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
        yield reaction_id

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
        line, _, comment = line.partition('#')
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
            if reaction.properties[prop] is not None:
                d[prop] = reaction.properties[prop]
            if prop == 'equation' and not is_equation_valid(d[prop]):
                del d[prop]

        return d

    def write_compounds(self, stream, compounds, properties=None):
        """Write iterable of compounds as YAML object to stream.

        Args:
            stream: File-like object.
            compounds: Iterable of compound entries.
            properties: Set of compound properties to output (or None to output
                all).
        """
        def iter_entries():
            for c in compounds:
                entry = self.convert_compound_entry(c)
                if entry is None:
                    continue
                if properties is not None:
                    entry = OrderedDict(
                        (key, value) for key, value in iteritems(entry)
                        if key == 'id' or key in properties)
                yield entry

        self._dump(stream, list(iter_entries()))

    def write_reactions(self, stream, reactions, properties=None):
        """Write iterable of reactions as YAML object to stream.

        Args:
            stream: File-like object.
            compounds: Iterable of reaction entries.
            properties: Set of reaction properties to output (or None to output
                all).
        """
        def iter_entries():
            for r in reactions:
                entry = self.convert_reaction_entry(r)
                if entry is None:
                    continue
                if properties is not None:
                    entry = OrderedDict(
                        (key, value) for key, value in iteritems(entry)
                        if key == 'id' or key in properties)
                yield entry

        self._dump(stream, list(iter_entries()))
