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

"""Module for reading and writing native formats.

These formats are either table-based or YAML-based. Table-based formats
are space-separated and empty lines are ignored. Comments starting with
pound (#). YAML-based formats are structured data following the YAML
specification.
"""

from __future__ import absolute_import

import os
import logging
import re
import csv

import yaml
from six import string_types, iteritems

from ..reaction import Reaction, Compound
from .context import FilePathContext, FileMark
from . import modelseed


# Module-level logging
logger = logging.getLogger(__name__)

# Model files to try to open if a directory was specified
DEFAULT_MODEL = ('model.yaml', 'model.yml')


class ParseError(Exception):
    """Exception used to signal errors while parsing"""


def whendefined(func, value):
    """Apply func to value if value is not None"""
    return func(value) if value is not None else None


class CompoundEntry(object):
    """Representation of a compound entry in a native model"""

    def __init__(self, compound_id, properties, filemark=None):
        self._id = compound_id
        self._properties = dict(properties)
        self._filemark = filemark

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        return self._properties.get('name')

    @property
    def formula(self):
        return self._properties.get('formula')

    @property
    def charge(self):
        return whendefined(int, self._properties.get('charge'))

    @property
    def kegg(self):
        return self._properties.get('kegg')

    @property
    def cas(self):
        return self._properties.get('cas')

    @property
    def zeromass(self):
        return self._properties.get('zeromass')

    @property
    def properties(self):
        return self._properties

    @property
    def filemark(self):
        return self._filemark


class ReactionEntry(object):
    """Representation of a reaction entry in a native model"""

    def __init__(self, id, properties, filemark=None):
        self._id = id
        self._properties = dict(properties)
        self._name = self._properties.get('name')
        self._equation = self._properties.get('equation')
        self._ec = self._properties.get('ec')
        self._genes = self._properties.get('genes')
        self._filemark = filemark

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        return self._name

    @property
    def equation(self):
        return self._equation

    @property
    def ec(self):
        return self._ec

    @property
    def genes(self):
        return self._genes

    @property
    def properties(self):
        return self._properties

    @property
    def filemark(self):
        return self._filemark


class NativeModel(object):
    """Represents a model specified using the native data formats

    The model is created from a model file or from a directory containing a
    model file using the default file name (model.yaml or model.yml). This file
    can specify the model fully or refer to other files within the same
    directory subtree that specifes part of the model.
    """

    def __init__(self, path):
        """Create a model from the specified model file or directory"""

        if os.path.isfile(path):
            self._context = FilePathContext(path)
            with open(self._context.filepath, 'r') as f:
                self._model = yaml.load(f)
        else:
            # Try to open the default file
            for filename in DEFAULT_MODEL:
                try:
                    self._context = FilePathContext(
                        os.path.join(path, filename))
                    with open(self._context.filepath, 'r') as f:
                        self._model = yaml.load(f)
                        break
                except Exception:
                    logger.debug('Failed to load model file', exc_info=True)
            else:
                # No model could be loaded
                raise ParseError('No model file could be found ({})'.format(
                    ', '.join(DEFAULT_MODEL)))

    def get_name(self):
        """Return the name specified by the model"""
        return self._model.get('name', None)

    def get_biomass_reaction(self):
        """Return the biomass reaction specified by the model"""
        return self._model.get('biomass', None)

    def get_default_flux_limit(self):
        """Return the default flux limit specified by the model"""
        return self._model.get('default_flux_limit', None)

    def parse_reactions(self):
        """Yield tuples of reaction ID and reactions defined in the model"""

        # Parse reactions defined in the main model file
        if 'reactions' in self._model:
            for reaction in parse_reaction_list(
                    self._context, self._model['reactions']):
                yield reaction

    def parse_model(self):
        """Yield reaction IDs of model reactions"""

        if 'model' in self._model:
            for reaction_id in parse_model_group_list(
                    self._context, self._model['model']):
                yield reaction_id
        else:
            reactions = set(reaction.id for reaction in self.parse_reactions())
            for reaction in reactions:
                yield reaction

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

        if 'media' in self._model:
            if not isinstance(self._model['media'], list):
                raise ParseError('Expected media to be a list')

            for medium_compound in parse_medium_list(
                    self._context, self._model['media']):
                yield medium_compound

    def parse_compounds(self):
        """Yield CompoundEntries for defined compounds"""

        if 'compounds' in self._model:
            for compound in parse_compound_list(
                    self._context, self._model['compounds']):
                yield compound

    @property
    def context(self):
        return self._context


def parse_compound(compound_def, context=None):
    """Parse a structured compound definition as obtained from a YAML file

    Returns a CompoundEntry."""

    compound_id = compound_def.get('id')
    if compound_id is None:
        raise ParseError('Compound ID missing')

    mark = FileMark(context, None, None)
    return CompoundEntry(compound_id, compound_def, mark)


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

    for i, row in enumerate(csv.DictReader(f, delimiter='\t')):
        if 'id' not in row or row['id'].strip() == '':
            raise ParseError('Expected `id` column in table')

        props = {key: value for key, value in iteritems(row) if value != ''}
        mark = FileMark(context, i + 2, None)
        yield CompoundEntry(row['id'], props, mark)


def parse_compound_yaml_file(path, f):
    """Parse a file as a YAML-format list of compounds

    Path can be given as a string or a context.
    """

    return parse_compound_list(path, yaml.load(f))


def parse_compound_file(path, format):
    """Open and parse reaction file based on file extension or given format

    Path can be given as a string or a context.
    """

    context = FilePathContext(path)

    # YAML files do not need to explicitly specify format
    if (re.match(r'.+\.(yml|yaml)$', context.filepath) and
            (format is None or format == 'yaml')):
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
    elif re.match(r'.+\.tsv$', context.filepath) or format == 'tsv':
        logger.debug('Parsing compound file {} as TSV'.format(
            context.filepath))
        with context.open('r') as f:
            for compound in parse_compound_table_file(context, f):
                yield compound
    else:
        raise ParseError('Unable to detect format of compound file {}'.format(
            context.filepath))


def parse_reaction_equation(equation_def):
    """Parse a structured reaction equation as obtained from a YAML file

    Returns a Reaction.
    """

    def parse_compound_list(l):
        """Parse a list of reactants or metabolites"""
        for compound_def in l:
            compound_id = compound_def.get('id')
            if compound_id is None:
                raise ParseError('Compound ID missing')

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
        return modelseed.parse_reaction(equation_def).normalized()

    compartment = equation_def.get('compartment', None)
    reversible = bool(equation_def.get('reversible', True))
    left = equation_def.get('left', [])
    right = equation_def.get('right', [])
    if len(left) == 0 and len(right) == 0:
        raise ParseError('Reaction values are missing')

    return Reaction(Reaction.Bidir if reversible else Reaction.Right,
                    parse_compound_list(left), parse_compound_list(right))


def parse_reaction(reaction_def, context=None):
    """Parse a structured reaction definition as obtained from a YAML file

    Returns a ReactionEntry.
    """

    reaction_id = reaction_def.get('id')
    if reaction_id is None:
        raise ParseError('Reaction ID missing')

    reaction_props = dict(reaction_def)

    # Parse reaction equation
    if 'equation' in reaction_def:
        reaction_props['equation'] = (
            parse_reaction_equation(reaction_def['equation']))

    mark = FileMark(context, None, None)
    return ReactionEntry(reaction_id, reaction_props, mark)


def parse_reaction_list(path, reactions):
    """Parse a structured list of reactions as obtained from a YAML file

    Yields tuples of reaction ID and reaction object. Path can be given as a
    string or a context.
    """

    context = FilePathContext(path)

    for reaction_def in reactions:
        if 'include' in reaction_def:
            include_context = context.resolve(reaction_def['include'])
            for reaction in parse_reaction_file(include_context):
                yield reaction
        else:
            yield parse_reaction(reaction_def, context)


def parse_reaction_yaml_file(path, f):
    """Parse a file as a YAML-format list of reactions

    Path can be given as a string or a context.
    """

    return parse_reaction_list(path, yaml.load(f))


def parse_reaction_table_file(path, f):
    """Parse a tab-separated file containing reaction IDs and properties

    The reaction properties are parsed according to the header which specifies
    which property is contained in each column.
    """

    context = FilePathContext(path)

    for lineno, row in enumerate(csv.DictReader(f, delimiter='\t')):
        if 'id' not in row or row['id'].strip() == '':
            raise ParseError('Expected `id` column in table')

        props = {key: value for key, value in iteritems(row) if value != ''}

        if 'equation' in props:
            props['equation'] = modelseed.parse_reaction(props['equation'])

        mark = FileMark(context, lineno + 2, 0)
        yield ReactionEntry(row['id'], props, mark)


def parse_reaction_file(path):
    """Open and parse reaction file based on file extension

    Path can be given as a string or a context.
    """

    context = FilePathContext(path)

    if re.match(r'.+\.tsv$', context.filepath):
        logger.debug('Parsing reaction file {} as TSV'.format(
            context.filepath))
        with context.open('r') as f:
            for reaction in parse_reaction_table_file(context, f):
                yield reaction
    elif re.match(r'.+\.(yml|yaml)$', context.filepath):
        logger.debug('Parsing reaction file {} as YAML'.format(
            context.filepath))
        with context.open('r') as f:
            for reaction in parse_reaction_yaml_file(context, f):
                yield reaction
    else:
        raise ParseError('Unable to detect format of reaction file {}'.format(
            context.filepath))


def parse_medium(medium_def):
    """Parse a structured medium definition as obtained from a YAML file

    Returns in iterator of compound, reaction, lower and upper bounds.
    """

    default_compartment = medium_def.get('compartment')

    for compound_def in medium_def.get('compounds', []):
        compartment = compound_def.get('compartment', default_compartment)
        compound = Compound(compound_def['id'], compartment=compartment)
        reaction = compound_def.get('reaction')
        lower = compound_def.get('lower')
        upper = compound_def.get('upper')
        yield compound, reaction, lower, upper


def parse_medium_list(path, medium):
    """Parse a structured medium list as obtained from a YAML file.

    Yields tuples of compound, reaction ID, lower and upper flux bounds. Path
    can be given as a string or a context.
    """

    context = FilePathContext(path)

    for medium_def in medium:
        if 'include' in medium_def:
            include_context = context.resolve(medium_def['include'])
            for medium_compound in parse_medium_file(include_context):
                yield medium_compound
        else:
            for medium_compound in parse_medium(medium_def):
                yield medium_compound


def parse_medium_yaml_file(path, f):
    """Parse a file as a YAML-format medium definition

    Path can be given as a string or a context.
    """

    return parse_medium(yaml.load(f))


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


def parse_medium_file(path):
    """Parse a file as a list of medium compounds with flux limits

    The file format is detected and the file is parsed accordingly. Path can
    be given as a string or a context.
    """

    context = FilePathContext(path)

    if re.match(r'.+\.tsv$', context.filepath):
        logger.debug('Parsing medium file {} as TSV'.format(
            context.filepath))
        with context.open('r') as f:
            for entry in parse_medium_table_file(f):
                yield entry
    elif re.match(r'.+\.(yml|yaml)$', context.filepath):
        logger.debug('Parsing medium file {} as YAML'.format(
            context.filepath))
        with context.open('r') as f:
            for entry in parse_medium_yaml_file(context, f):
                yield entry
    else:
        raise ParseError('Unable to detect format of medium file {}'.format(
            context.filepath))


def parse_limit(limit_def):
    """Parse a structured flux limit definition as obtained from a YAML file

    Returns a tuple of reaction, lower and upper bound.
    """

    if 'reaction' not in limit_def:
        raise ParseError('Expected reaction key in limit entry')

    reaction = limit_def.get('reaction')
    lower = limit_def.get('lower', None)
    upper = limit_def.get('upper', None)

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

    return parse_limits_list(path, yaml.load(f))


def parse_limits_file(path):
    """Parse a file as a list of reaction flux limits

    The file format is detected and the file is parsed accordingly. Path can
    be given as a string or a context.
    """

    context = FilePathContext(path)

    if re.match(r'.+\.tsv$', context.filepath):
        logger.debug('Parsing limits file {} as TSV'.format(
            context.filepath))
        with context.open('r') as f:
            for limit in parse_limits_table_file(f):
                yield limit
    elif re.match(r'.+\.(yml|yaml)$', context.filepath):
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
    return parse_model_group_list(path, yaml.load(f))


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

    if re.match(r'.+\.tsv$', context.filepath):
        logger.debug('Parsing model file {} as TSV'.format(context.filepath))
        with context.open('r') as f:
            for reaction_id in parse_model_table_file(context, f):
                yield reaction_id
    elif re.match(r'.+\.(yml|yaml)$', context.filepath):
        logger.debug('Parsing model file {} as YAML'.format(context.filepath))
        with context.open('r') as f:
            for reaction_id in parse_model_yaml_file(context, f):
                yield reaction_id
