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

"""Parser for SBML model files."""

from __future__ import unicode_literals

from decimal import Decimal
from fractions import Fraction
from functools import partial
import logging
import re
import json
from collections import OrderedDict, Counter

# Import ElementTree XML parser. The lxml etree implementation may also be
# used with SBMLReader but has compatibility issues with SBMLWriter.
try:
    import xml.etree.cElementTree as ETree
except ImportError:
    import xml.etree.ElementTree as ETree

from six import itervalues, iteritems, text_type, PY3

from .context import FileMark
from .entry import (CompoundEntry as BaseCompoundEntry,
                    ReactionEntry as BaseReactionEntry,
                    CompartmentEntry as BaseCompartmentEntry,
                    DictReactionEntry, DictCompoundEntry,
                    DictCompartmentEntry)
from .native import NativeModel
from ..reaction import Reaction, Compound, Direction
from ..metabolicmodel import create_exchange_id
from ..expression.boolean import Expression, And, Or, Variable
from .. import util

logger = logging.getLogger(__name__)

# Level 1 namespaces
SBML_NS_L1 = 'http://www.sbml.org/sbml/level1'

# Level 2 namespaces
SBML_NS_L2 = 'http://www.sbml.org/sbml/level2'
SBML_NS_L2_V2 = 'http://www.sbml.org/sbml/level2/version2'
SBML_NS_L2_V3 = 'http://www.sbml.org/sbml/level2/version3'
SBML_NS_L2_V4 = 'http://www.sbml.org/sbml/level2/version4'
SBML_NS_L2_V5 = 'http://www.sbml.org/sbml/level2/version5'

# Level 3 namespaces
SBML_NS_L3_V1_CORE = 'http://www.sbml.org/sbml/level3/version1/core'

MATHML_NS = 'http://www.w3.org/1998/Math/MathML'
XHTML_NS = 'http://www.w3.org/1999/xhtml'

# FBC namespaces
FBC_V1 = 'http://www.sbml.org/sbml/level3/version1/fbc/version1'
FBC_V2 = 'http://www.sbml.org/sbml/level3/version1/fbc/version2'


# COBRA ID mappings
_COBRA_DECODE_ESCAPES = {
    '_DASH_': '-',
    '_FSLASH_': '/',
    '_BSLASH_': '\\',
    '_LPAREN_': '(',
    '_RPAREN_': ')',
    '_LSQBKT_': '[',
    '_RSQBKT_': ']',
    '_COMMA_': ',',
    '_PERIOD_': '.',
    '_APOS_': "'"
}


def _tag(tag, namespace=None):
    """Prepend namespace to tag name"""
    if namespace is None:
        return text_type(tag)
    return '{{{}}}{}'.format(namespace, tag)


class ParseError(Exception):
    """Error parsing SBML file"""


class _SBMLEntry(object):
    """Base class for compound and reaction entries."""

    def __init__(self, reader, root):
        self._reader = reader
        self._root = root
        self._id = self._element_get_id(root)

    def _element_get_id(self, element):
        """Get id of reaction or species element.

        In old levels the name is used as the id. This method returns the
        correct attribute depending on the level.
        """
        if self._reader._level > 1:
            entry_id = element.get('id')
        else:
            entry_id = element.get('name')
        return entry_id

    @property
    def id(self):
        """Entity ID"""
        return self._id

    @property
    def xml_notes(self):
        """Access the entity notes as an XML document fragment."""
        return self._root.find(self._reader._sbml_tag('notes'))

    @property
    def xml_annotation(self):
        """Access the entity annotation as an XML document fragment."""
        return self._root.find(self._reader._sbml_tag('annotation'))


class SBMLSpeciesEntry(_SBMLEntry, BaseCompoundEntry):
    """Species entry in the SBML file"""

    def __init__(self, reader, root, filemark=None):
        super(SBMLSpeciesEntry, self).__init__(reader, root)

        self._name = root.get('name')
        self._comp = root.get('compartment')

        if self._comp is None:
            msg = 'Species {} has no compartment!'.format(self.id)
            if self._reader._strict:
                raise ParseError(msg)
            else:
                logger.warning(msg)

        self._boundary = root.get('boundaryCondition', 'false') == 'true'

        # In non-strict mode the species that ends with _b are considered
        # boundary conditions.
        if not self._reader._strict and self._id.endswith('_b'):
            logger.warning('Species {} was converted to boundary condition'
                           ' because of "_b" suffix'.format(self.id))
            self._boundary = True

        self._filemark = filemark

    @property
    def name(self):
        """Species name"""
        return self._name

    @property
    def compartment(self):
        """Species compartment"""
        return self._comp

    def _parse_charge_string(self, s):
        try:
            return int(s)
        except ValueError:
            if self._reader._strict:
                raise ParseError('Invalid charge for species {}'.format(
                    self.id))
            else:
                return None

    @property
    def charge(self):
        """Species charge"""
        if self._reader._level == 3:
            # Look for FBC charge
            for ns in (FBC_V2, FBC_V1):
                charge = self._root.get(_tag('charge', ns))
                if charge is not None:
                    return self._parse_charge_string(charge)
        else:
            charge = self._root.get('charge')
            if charge is not None:
                return self._parse_charge_string(charge)

        return None

    @property
    def formula(self):
        """Species formula"""
        if self._reader._level == 3:
            for ns in (FBC_V2, FBC_V1):
                formula = self._root.get(_tag('chemicalFormula', ns))
                if formula is not None:
                    return formula

        return None

    @property
    def boundary(self):
        """Whether this compound is a boundary condition"""
        return self._boundary

    @property
    def properties(self):
        """All species properties as a dict"""
        properties = {'id': self._id,
                      'boundary': self._boundary}
        if 'name' in self._root.attrib:
            properties['name'] = self._root.get('name')
        if 'compartment' in self._root.attrib:
            properties['compartment'] = self._root.get('compartment')

        charge = self.charge
        if charge is not None:
            properties['charge'] = charge

        formula = self.formula
        if formula is not None:
            properties['formula'] = formula

        return properties

    @property
    def filemark(self):
        return self._filemark


class SBMLReactionEntry(_SBMLEntry, BaseReactionEntry):
    """Reaction entry in SBML file"""

    def __init__(self, reader, root, filemark=None):
        super(SBMLReactionEntry, self).__init__(reader, root)

        self._name = self._root.get('name')
        self._rev = self._root.get('reversible', 'true') == 'true'

        left, right = [], []
        for side, tag_name in ((left, 'listOfReactants'),
                               (right, 'listOfProducts')):
            for species_id, value in self._parse_species_references(tag_name):
                try:
                    species_entry = self._reader.get_species(species_id)
                    if (self._reader._ignore_boundary and
                            species_entry.boundary):
                        # Skip boundary species when ignoring these
                        continue
                except KeyError:
                    message = ('Reaction {} references non-existent'
                               ' species {}'.format(self._id, species_id))
                    if not self._reader._strict:
                        # In non-strict mode simply skip these references
                        logger.warn(message)
                        continue
                    raise ParseError(message)

                species_id, species_comp = (
                    species_entry.id, species_entry.compartment)
                compound = Compound(species_id, compartment=species_comp)
                side.append((compound, value))

        direction = Direction.Both if self._rev else Direction.Forward
        self._equation = Reaction(direction, left, right)

        # Parse flux bounds of reaction
        self._lower_flux = None
        self._upper_flux = None

        if reader._level == 3:
            lower = self._root.get(_tag('lowerFluxBound', FBC_V2))
            if lower is not None:
                if lower not in reader._model_constants:
                    raise ParseError(
                        'Lower flux parameter of {} is not defined in'
                        ' model parameter constants'.format(self.id))
                self._lower_flux = reader._model_constants[lower]

            upper = self._root.get(_tag('upperFluxBound', FBC_V2))
            if upper is not None:
                if upper not in reader._model_constants:
                    raise ParseError(
                        'Upper flux parameter of {} is not defined in'
                        ' model parameter constants'.format(self.id))
                self._upper_flux = reader._model_constants[upper]

            if (lower is None and upper is None and
                    self.id in reader._reaction_flux_bounds):
                # Parse bounds from listOfFluxBounds in FBCv1
                for bound in reader._reaction_flux_bounds[self.id]:
                    if bound.operation == 'equal':
                        self._lower_flux = bound.value
                        self._upper_flux = bound.value
                    elif bound.operation == 'lessEqual':
                        self._upper_flux = bound.value
                    elif bound.operation == 'greaterEqual':
                        self._lower_flux = bound.value

        self._filemark = filemark

    def _parse_species_references(self, name):
        """Yield species id and parsed value for a speciesReference list"""
        for species in self._root.iterfind('./{}/{}'.format(
                self._reader._sbml_tag(name),
                self._reader._sbml_tag('speciesReference'))):

            species_id = species.get('species')

            if self._reader._level == 1:
                # In SBML level 1 only positive integers are allowed for
                # stoichiometry but a positive integer denominator can be
                # specified.
                try:
                    value = int(species.get('stoichiometry', 1))
                    denom = int(species.get('denominator', 1))
                    species_value = Fraction(value, denom)
                except ValueError:
                    message = ('Non-integer stoichiometry is not allowed in'
                               ' SBML level 1 (reaction {})'.format(self.id))
                    if self._reader._strict:
                        raise ParseError(message)
                    else:
                        logger.warning(message)
                    species_value = Decimal(species.get('stoichiometry', 1))
            elif self._reader._level == 2:
                # Stoichiometric value is a double but can alternatively be
                # specified using math (not implemented).
                value_str = species.get('stoichiometry', None)
                if value_str is None:
                    if 'stoichiometryMath' in species:
                        raise ParseError('stoichiometryMath in '
                                         'speciesReference is not implemented')
                    species_value = 1
                else:
                    species_value = Decimal(value_str)
            elif self._reader._level == 3:
                # Stoichiometric value is a double but can alternatively be
                # specified using initial assignment (not implemented).
                value_str = species.get('stoichiometry', None)
                if value_str is None:
                    raise ParseError('Missing stoichiometry in'
                                     ' speciesReference is not implemented')
                species_value = Decimal(value_str)

            if species_value % 1 == 0:
                species_value = int(species_value)

            yield species_id, species_value

    @property
    def id(self):
        """Reaction ID"""
        return self._id

    @property
    def name(self):
        """Reaction name"""
        return self._name

    @property
    def reversible(self):
        """Whether the reaction is reversible"""
        return self._rev

    @property
    def equation(self):
        """Reaction equation is a :class:`Reaction <psamm.reaction.Reaction>`
        object"""
        return self._equation

    @property
    def kinetic_law_reaction_parameters(self):
        """Iterator over the values of kinetic law reaction parameters"""

        for parameter in self._root.iterfind(
                './{}/{}/{}'.format(self._reader._sbml_tag('kineticLaw'),
                                    self._reader._sbml_tag('listOfParameters'),
                                    self._reader._sbml_tag('parameter'))):
            param_id = parameter.get('id')
            param_name = parameter.get('name')
            param_value = Decimal(parameter.get('value'))
            param_units = parameter.get('units')

            yield param_id, param_name, param_value, param_units

    @property
    def properties(self):
        """All reaction properties as a dict"""
        properties = {'id': self._id,
                      'reversible': self._rev,
                      'equation': self._equation}
        if 'name' in self._root.attrib:
            properties['name'] = self._root.get('name')
        if self._lower_flux is not None:
            properties['lower_flux'] = self._lower_flux
        if self._upper_flux is not None:
            properties['upper_flux'] = self._upper_flux

        return properties

    @property
    def filemark(self):
        return self._filemark


class SBMLCompartmentEntry(_SBMLEntry, BaseCompartmentEntry):
    """Compartment entry in the SBML file"""

    def __init__(self, reader, root, filemark=None):
        super(SBMLCompartmentEntry, self).__init__(reader, root)

        self._name = self._root.get('name')
        self._filemark = filemark

    @property
    def properties(self):
        """All compartment properties as a dict."""
        properties = {'id': self._id}
        if self._name is not None:
            properties['name'] = self._name

        return properties

    @property
    def filemark(self):
        return self._filemark


class SBMLObjectiveEntry(object):
    """Flux objective defined with FBC"""

    def __init__(self, reader, namespace, root):
        self._reader = reader
        self._root = root
        self._namespace = namespace

        self._id = self._root.get(_tag('id', self._namespace))
        if self._id is None:
            raise ParseError('Objective has no "id" attribute')

        self._name = self._root.get(_tag('name', self._namespace))
        self._type = self._root.get(_tag('type', self._namespace))

        if self._type is None:
            raise ParseError(
                'Missing "type" attribute on objective: {}'.format(self.id))

        # Find flux objectives
        self._reactions = {}
        for fo in self._root.iterfind('./{}/{}'.format(
                _tag('listOfFluxObjectives', self._namespace),
                _tag('fluxObjective', self._namespace))):
            reaction = fo.get(_tag('reaction', self._namespace))
            coefficient = Decimal(fo.get(_tag('coefficient', self._namespace)))
            if reaction is None or coefficient is None:
                raise ParseError(
                    'Missing attributes on flux objective: {}'.format(self.id))

            if coefficient != 0:
                self._reactions[reaction] = coefficient

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        return self._name

    @property
    def type(self):
        return self._type

    @property
    def reactions(self):
        return iteritems(self._reactions)


class SBMLFluxBoundEntry(object):
    """Flux bound defined with FBC V1.

    Flux bounds defined with FBC V2 are instead encoded as ``upper_flux`` and
    ``lower_flux`` properties on the ReactionEntry objects.
    """

    LESS_EQUAL = 'lessEqual'
    GREATER_EQUAL = 'greaterEqual'
    EQUAL = 'equal'

    def __init__(self, reader, namespace, root):
        self._reader = reader
        self._root = root
        self._namespace = namespace

        self._id = self._root.get(_tag('id', self._namespace))
        self._name = self._root.get(_tag('name', self._namespace))

        self._reaction = self._root.get(_tag('reaction', self._namespace))
        if self._reaction is None:
            raise ParseError('Flux bound is missing "reaction" attribute')

        self._operation = self._root.get(_tag('operation', self._namespace))
        if self._operation is None:
            raise ParseError('Flux bound is missing "operation" attribute')
        elif self._operation not in {
                self.LESS_EQUAL, self.GREATER_EQUAL, self.EQUAL}:
            raise ParseError('Invalid flux bound operation: {}'.format(
                self._operation))

        value = self._root.get(_tag('value', self._namespace))
        if value is None:
            raise ParseError('Flux bound is missing "value" attribute')
        try:
            self._value = Decimal(value)
        except ValueError:
            raise ParseError('Unable to parse flux bound value')

    @property
    def id(self):
        """Return ID of flux bound."""
        return self._id

    @property
    def name(self):
        """Return name of flux bound."""
        return self._name

    @property
    def reaction(self):
        """Return reaction ID that the flux bound pertains to."""
        return self._reaction

    @property
    def operation(self):
        """Return the operation of the flux bound.

        Returns one of LESS_EQUAL, GREATER_EQUAL or EQUAL.
        """
        return self._operation

    @property
    def value(self):
        """Return the flux bound value."""
        return self._value


class SBMLReader(object):
    """Reader of SBML model files

    The constructor takes a file-like object which will be parsed as XML and
    then as SBML according to the specification. If the ``strict`` parameter is
    set to False, the parser will revert to a more lenient parsing which is
    required for many older models. This tries to mimic the inconsistencies
    employed by COBRA when parsing models.

    If ``ignore_boundary`` is ``True``, the species that are marked as
    boundary conditions will simply be dropped from the species list and from
    the reaction equations, and any boundary compartment will be dropped too.
    Otherwise the boundary species will be retained. Retaining these is only
    useful to extract specific information from those species objects.

    Args:
        file: File-like object to parse XML SBML content from.
        strict: Indicating whether strict parsing is enabled.
        ignore_boundary: Indicating whether boundary species are dropped.
        context: Optional file parsing context from
            :mod:`psamm.datasource.context`.
    """

    def __init__(self, file, strict=False, ignore_boundary=True,
                 context=None):
        # Parse SBML file
        tree = ETree.parse(file)
        root = tree.getroot()

        self._strict = strict
        self._ignore_boundary = ignore_boundary
        self._context = context

        # Parse level and version
        self._sbml_tag = None
        self._level = int(root.get('level'))
        self._version = int(root.get('version'))

        if self._level == 1:
            self._sbml_tag = partial(_tag, namespace=SBML_NS_L1)
        elif self._level == 2:
            if self._version == 1:
                self._sbml_tag = partial(_tag, namespace=SBML_NS_L2)
            elif self._version == 2:
                self._sbml_tag = partial(_tag, namespace=SBML_NS_L2_V2)
            elif self._version == 3:
                self._sbml_tag = partial(_tag, namespace=SBML_NS_L2_V3)
            elif self._version == 4:
                self._sbml_tag = partial(_tag, namespace=SBML_NS_L2_V4)
            elif self._version == 5:
                self._sbml_tag = partial(_tag, namespace=SBML_NS_L2_V5)
            else:
                if self._strict:
                    raise ParseError('SBML level 2, version {}'
                                     ' not implemented'.format(self._version))
                # Assume latest version
                self._version = 5
                self._sbml_tag = partial(_tag, namespace=SBML_NS_L2_V5)
        elif self._level == 3:
            if self._version == 1:
                self._sbml_tag = partial(_tag, namespace=SBML_NS_L3_V1_CORE)
            else:
                if self._strict:
                    raise ParseError('SBML level 3, version {}'
                                     ' not implemented'.format(self._version))
                # Assume latest version
                self._version = 1
                self._sbml_tag = partial(_tag, namespace=SBML_NS_L3_V1_CORE)
        else:
            if self._strict:
                raise ParseError('SBML level {} not implemented'.format(
                    self._level))
            # Assume level 1
            self._level = 1
            self._sbml_tag = partial(_tag, namespace=SBML_NS_L1)

        self._model = root.find(self._sbml_tag('model'))

        # Parameters
        self._model_constants = {}
        params = self._model.find(self._sbml_tag('listOfParameters'))
        if params is not None:
            for param in params.iterfind(self._sbml_tag('parameter')):
                if param.get('constant') == 'true':
                    param_id = param.get('id')
                    value = Decimal(param.get('value'))
                    self._model_constants[param_id] = value

        # Flux bounds
        self._flux_bounds = []
        self._reaction_flux_bounds = {}
        if self._level == 3:
            # Only represented with this tag in FBC V1. In FBC V2 the flux
            # bounds are instead represented in the global listOfParameters
            # and referenced by attributes on the reactions. Only with FBC V1
            # are FluxBoundEntry objects created.
            flux_bounds = self._model.find(_tag('listOfFluxBounds', FBC_V1))
            if flux_bounds is not None:
                for flux_bound in flux_bounds.iterfind(
                        _tag('fluxBound', FBC_V1)):
                    entry = SBMLFluxBoundEntry(self, FBC_V1, flux_bound)
                    self._flux_bounds.append(entry)

                    # Create reference from reaction to flux bound
                    entries = self._reaction_flux_bounds.setdefault(
                        entry.reaction, [])
                    entries.append(entry)

        # Compartments
        self._model_compartments = {}
        self._compartments = self._model.find(
            self._sbml_tag('listOfCompartments'))
        for compartment in self._compartments.iterfind(
                self._sbml_tag('compartment')):
            filemark = FileMark(
                self._context, self._get_sourceline(compartment), None)
            entry = SBMLCompartmentEntry(self, compartment, filemark=filemark)
            self._model_compartments[entry.id] = entry

        # Species
        self._model_species = {}
        self._species = self._model.find(self._sbml_tag('listOfSpecies'))
        for species in self._species.iterfind(self._sbml_tag('species')):
            filemark = FileMark(
                self._context, self._get_sourceline(species), None)
            entry = SBMLSpeciesEntry(self, species, filemark=filemark)
            self._model_species[entry.id] = entry

        # Reactions
        self._model_reactions = {}
        self._reactions = self._model.find(self._sbml_tag('listOfReactions'))
        for reaction in self._reactions.iterfind(self._sbml_tag('reaction')):
            filemark = FileMark(
                self._context, self._get_sourceline(reaction), None)
            entry = SBMLReactionEntry(self, reaction, filemark=filemark)
            self._model_reactions[entry.id] = entry

        if self._ignore_boundary:
            # Remove compartments that only contain boundary compounds
            empty_compartments = set(self._model_compartments)
            valid_compartments = set()
            for species in itervalues(self._model_species):
                empty_compartments.discard(species.compartment)
                if not species.boundary:
                    valid_compartments.add(species.compartment)

            boundary_compartments = (
                set(self._model_compartments) - empty_compartments -
                valid_compartments)
            for compartment in boundary_compartments:
                logger.info('Ignoring boundary compartment: {}'.format(
                    compartment))
                del self._model_compartments[compartment]

        # Objectives
        self._model_objectives = {}
        self._active_objective = None
        if self._level == 3:
            objectives = None
            for fbc_ns in (FBC_V2, FBC_V1):
                objectives = self._model.find(_tag('listOfObjectives', fbc_ns))
                if objectives is not None:
                    break

            if objectives is not None:
                for objective in objectives.iterfind(
                        _tag('objective', fbc_ns)):
                    entry = SBMLObjectiveEntry(self, fbc_ns, objective)
                    self._model_objectives[entry.id] = entry

                active = objectives.get(_tag('activeObjective', fbc_ns))
                if active is None or active not in self._model_objectives:
                    raise ParseError('Active objective is invalid: {}'.format(
                        active))

                self._active_objective = self._model_objectives[active]

    def _get_sourceline(self, element):
        """Get source line of element (only supported by lxml)."""
        return getattr(element, 'sourceline', None)

    def get_compartment(self, compartment_id):
        """Return :class:`.CompartmentEntry` corresponding to id."""
        return self._model_compartments[compartment_id]

    def get_reaction(self, reaction_id):
        """Return :class:`.SBMLReactionEntry` corresponding to reaction_id"""
        return self._model_reactions[reaction_id]

    def get_species(self, species_id):
        """Return :class:`.SBMLSpeciesEntry` corresponding to species_id"""
        return self._model_species[species_id]

    def get_objective(self, objective_id):
        """Return :class:`.SBMLObjectiveEntry` corresponding to objective_id"""
        return self._model_objectives[objective_id]

    def get_active_objective(self):
        return self._active_objective

    @property
    def compartments(self):
        """Iterator over :class:`.SBMLCompartmentEntry` entries."""
        return itervalues(self._model_compartments)

    @property
    def reactions(self):
        """Iterator over :class:`ReactionEntries <.SBMLReactionEntry>`"""
        return itervalues(self._model_reactions)

    @property
    def species(self):
        """Iterator over :class:`SpeciesEntries <.SpeciesEntry>`

        This will not yield boundary condition species if those are ignored.
        """
        return (c for c in itervalues(self._model_species)
                if not self._ignore_boundary or not c.boundary)

    @property
    def objectives(self):
        """Iterator over :class:`.SBMLObjectiveEntry`"""
        return itervalues(self._model_objectives)

    @property
    def flux_bounds(self):
        """Iterator over :class:`.SBMLFluxBoundEntry`"""
        return iter(self._flux_bounds)

    @property
    def id(self):
        """Model ID"""
        return self._model.get('id', None)

    @property
    def name(self):
        """Model name"""
        return self._model.get('name', None)

    def create_model(self):
        """Create model from reader.

        Returns:
            :class:`psamm.datasource.native.NativeModel`.
        """
        properties = {
            'name': self.name,
            'default_flux_limit': 1000
        }

        # Load objective as biomass reaction
        objective = self.get_active_objective()
        if objective is not None:
            reactions = dict(objective.reactions)
            if len(reactions) == 1:
                reaction, value = next(iteritems(reactions))
                if ((value < 0 and objective.type == 'minimize') or
                        (value > 0 and objective.type == 'maximize')):
                    properties['biomass'] = reaction

        model = NativeModel(properties)

        # Load compartments into model
        for compartment in self.compartments:
            model.compartments.add_entry(compartment)

        # Load compounds into model
        for compound in self.species:
            model.compounds.add_entry(compound)

        # Load reactions into model
        for reaction in self.reactions:
            model.reactions.add_entry(reaction)

        # Create model reaction set
        for reaction in model.reactions:
            model.model[reaction.id] = None

        # Convert reaction limits properties to proper limits
        for reaction in model.reactions:
            props = reaction.properties
            if 'lower_flux' in props or 'upper_flux' in props:
                lower = props.get('lower_flux')
                upper = props.get('upper_flux')
                model.limits[reaction.id] = reaction.id, lower, upper

        # Load model limits from FBC V1 bounds if present, i.e. if FBC V1 is
        # used instead of V2.
        limits_lower = {}
        limits_upper = {}
        for bounds in self.flux_bounds:
            reaction = bounds.reaction
            if reaction in model.limits:
                continue

            if bounds.operation == SBMLFluxBoundEntry.LESS_EQUAL:
                if reaction not in limits_upper:
                    limits_upper[reaction] = bounds.value
                else:
                    raise ParseError(
                        'Conflicting bounds for {}'.format(reaction))
            elif bounds.operation == SBMLFluxBoundEntry.GREATER_EQUAL:
                if reaction not in limits_lower:
                    limits_lower[reaction] = bounds.value
                else:
                    raise ParseError(
                        'Conflicting bounds for {}'.format(reaction))
            elif bounds.operation == SBMLFluxBoundEntry.EQUAL:
                if (reaction not in limits_lower and
                        reaction not in limits_upper):
                    limits_lower[reaction] = bounds.value
                    limits_upper[reaction] = bounds.value
                else:
                    raise ParseError(
                        'Conflicting bounds for {}'.format(reaction))

        for reaction in model.reactions:
            if reaction.id in limits_lower or reaction.id in limits_upper:
                lower = limits_lower.get(reaction.id, None)
                upper = limits_upper.get(reaction.id, None)
                model.limits[reaction.id] = reaction.id, lower, upper

        return model


class SBMLWriter(object):
    """Writer of SBML files."""

    def __init__(self, cobra_flux_bounds=False):
        self._namespace = SBML_NS_L3_V1_CORE
        self._sbml_tag = partial(_tag, namespace=self._namespace)
        self._cobra_flux_bounds = False

    def _make_safe_id(self, id):
        """Returns a modified id that has been made safe for SBML.

        Replaces or deletes the ones that aren't allowed.
        """

        substitutions = {
            '-': '_DASH_',
            '/': '_FSLASH_',
            '\\': '_BSLASH_',
            '(': '_LPAREN_',
            ')': '_RPAREN_',
            '[': '_LSQBKT_',
            ']': '_RSQBKT_',
            ',': '_COMMA_',
            '.': '_PERIOD_',
            "'": '_APOS_'
        }

        id = id.encode('ascii', 'backslashreplace').decode('ascii')
        id = re.sub(r'\(([a-z])\)$', '_\\1', id)
        for symbol, escape in iteritems(substitutions):
            id = id.replace(symbol, escape)
        id = re.sub(r'[^a-zA-Z0-9_]', '', id)
        return id

    def _get_flux_bounds(self, r_id, model, flux_limits, equation):
        """Read reaction's limits to set up strings for limits in the output file.
        """
        if r_id not in flux_limits or flux_limits[r_id][0] is None:
            if equation.direction == Direction.Forward:
                lower = 0
            else:
                lower = -model.default_flux_limit
        else:
            lower = flux_limits[r_id][0]

        if r_id not in flux_limits or flux_limits[r_id][1] is None:
            if equation.direction == Direction.Reverse:
                upper = 0
            else:
                upper = model.default_flux_limit
        else:
            upper = flux_limits[r_id][1]

        if lower % 1 == 0:
            lower = int(lower)
        if upper % 1 == 0:
            upper = int(upper)
        return text_type(lower), text_type(upper)

    def _make_safe_numerical_id(self, stri):
        return stri.replace('-', 'neg').replace('.', '_')

    def _add_gene_associations(self, r_id, r_genes, gene_ids, r_tag):
        """Adds all the different kinds of genes into a list."""
        genes = ETree.SubElement(
            r_tag, _tag('geneProductAssociation', FBC_V2))
        if isinstance(r_genes, list):
            e = Expression(And(*(Variable(i) for i in r_genes)))
        else:
            e = Expression(r_genes)
        gene_stack = [(e.root, genes)]
        while len(gene_stack) > 0:
            current, parent = gene_stack.pop()
            if isinstance(current, Or):
                gene_tag = ETree.SubElement(parent, _tag('or', FBC_V2))
            elif isinstance(current, And):
                gene_tag = ETree.SubElement(parent, _tag('and', FBC_V2))
            elif isinstance(current, Variable):
                gene_tag = ETree.SubElement(parent, _tag(
                    'geneProductRef', FBC_V2))
                if current.symbol not in gene_ids:
                    id = 'g_' + util.create_unique_id(
                        self._make_safe_id(current.symbol), gene_ids)
                    gene_ids[id] = current.symbol
                    gene_tag.set(_tag('geneProduct', FBC_V2), id)
            if isinstance(current, (Or, And)):
                for item in current:
                    gene_stack.append((item, gene_tag))

    def _add_fbc_objective(self, model_tag, obj_id):
        """Adds the objective(s) to the sbml document."""
        objective_list = ETree.SubElement(model_tag, _tag(
            'listOfObjectives', FBC_V2))
        objective_list.set(_tag('activeObjective', FBC_V2), 'O_1')
        objective_tag = ETree.SubElement(
            objective_list, _tag('objective', FBC_V2))
        objective_tag.set(_tag('id', FBC_V2), 'O_1')
        objective_tag.set(_tag('type', FBC_V2), 'maximize')
        flux_objective_list = ETree.SubElement(objective_tag, _tag(
            'listOfFluxObjectives', FBC_V2))
        flux_objective_tag = ETree.SubElement(flux_objective_list, _tag(
            'fluxObjective', FBC_V2))
        flux_objective_tag.set(_tag('reaction', FBC_V2), 'R_' + obj_id)
        flux_objective_tag.set(_tag('coefficient', FBC_V2), '1')

    def _add_gene_list(self, parent_tag, gene_id_dict):
        """Create list of all gene products as sbml readable elements."""
        list_all_genes = ETree.SubElement(parent_tag, _tag(
            'listOfGeneProducts', FBC_V2))
        for id, label in sorted(iteritems(gene_id_dict)):
            gene_tag = ETree.SubElement(
                list_all_genes, _tag('geneProduct', FBC_V2))
            gene_tag.set(_tag('id', FBC_V2), id)
            gene_tag.set(_tag('label', FBC_V2), label)

    def _add_properties_notes(self, parent_tag, properties):
        for prop, value in sorted(iteritems(properties)):
            p_tag = ETree.SubElement(parent_tag, _tag('p', XHTML_NS))
            try:
                s = json.dumps(value)
            except TypeError:
                s = json.dumps(text_type(value))
            p_tag.text = '{}: {}'.format(prop, s)

    def _indent(self, elem, level=0):
        i = "\n" + level*"  "
        if len(elem):
            if not elem.text or not elem.text.strip():
                elem.text = i + "  "
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
            for elem in elem:
                self._indent(elem, level+1)
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
        else:
            if level and (not elem.tail or not elem.tail.strip()):
                elem.tail = i

    def write_model(self, file, model, pretty=False):
        """Write a given model to file.

        Args:
            file: File-like object open for writing.
            model: Instance of :class:`NativeModel` to write.
            pretty: Whether to format the XML output for readability.
        """
        ETree.register_namespace('mathml', MATHML_NS)
        ETree.register_namespace('xhtml', XHTML_NS)
        ETree.register_namespace('fbc', FBC_V2)

        # Load compound information
        compound_name = {}
        compound_properties = {}
        for compound in model.compounds:
            compound_name[compound.id] = (
                compound.name if compound.name is not None else compound.id)
            compound_properties[compound.id] = compound.properties

        model_reactions = set(model.model)

        reaction_properties = {}
        biomass_id = None
        for r in model.reactions:
            if (model_reactions is not None and
                    r.id not in model_reactions):
                continue

            reaction_id = util.create_unique_id(
                self._make_safe_id(r.id), reaction_properties)
            if r.id == model.biomass_reaction:
                biomass_id = reaction_id

            reaction_properties[reaction_id] = r.properties

        # Add exchange reactions to reaction_properties,
        # also add flux limit info to flux_limits
        flux_limits = {}
        for compound, reaction_id, lower, upper in itervalues(model.exchange):
            # Create exchange reaction
            if reaction_id is None:
                reaction_id = create_exchange_id(reaction_properties, compound)
            reaction_id = util.create_unique_id(
                self._make_safe_id(reaction_id), reaction_properties)

            reaction_properties[reaction_id] = {
                'id': reaction_id,
                'equation': Reaction(Direction.Both, {compound: -1})
            }

            if lower is None:
                lower = -model.default_flux_limit
            if upper is None:
                upper = model.default_flux_limit
            flux_limits[reaction_id] = (lower, upper)

            # Create a dummy properties dict for undefined compounds
            if compound.name not in compound_properties:
                compound_properties[compound.name] = {
                    'id': compound.name
                }

        root = ETree.Element(self._sbml_tag('sbml'))
        root.set(self._sbml_tag('level'), '3')
        root.set(self._sbml_tag('version'), '1')
        root.set(_tag('required', FBC_V2), 'false')
        if model.version_string is not None:
            notes_tag = ETree.SubElement(root, self._sbml_tag('notes'))
            body_tag = ETree.SubElement(notes_tag, _tag('body', XHTML_NS))
            self._add_properties_notes(
                body_tag, {'model version': model.version_string})

        model_tag = ETree.SubElement(root, self._sbml_tag('model'))
        model_tag.set(_tag('strict', FBC_V2), 'true')
        if model.name is not None:
            model_tag.set(self._sbml_tag('name'), model.name)

        # Build mapping from Compound to species ID
        model_compartments = {}
        model_species = {}
        species_ids = set()
        for _, properties in iteritems(reaction_properties):
            for compound, _ in properties['equation'].compounds:
                if compound in model_species:
                    continue

                # Create a dummy properties dict for undefined compounds
                if compound.name not in compound_properties:
                    compound_properties[compound.name] = {
                        'id': compound.name
                    }

                compound_id = util.create_unique_id(
                    self._make_safe_id(compound.name), species_ids)
                model_species[compound] = compound_id
                species_ids.add(compound_id)
                if compound.compartment not in model_compartments:
                    model_compartments[
                        compound.compartment] = 'C_' + util.create_unique_id(
                            self._make_safe_id(compound.compartment),
                            model_compartments)

        # Create list of compartments
        compartments = ETree.SubElement(
            model_tag, self._sbml_tag('listOfCompartments'))
        for _, compartment_id in iteritems(model_compartments):
            compartment_tag = ETree.SubElement(
                compartments, self._sbml_tag('compartment'))
            compartment_tag.set(self._sbml_tag('id'), compartment_id)
            compartment_tag.set(self._sbml_tag('constant'), 'true')

        # Create list of species
        species_list = ETree.SubElement(
            model_tag, self._sbml_tag('listOfSpecies'))
        for species, species_id in sorted(
                iteritems(model_species), key=lambda x: x[1]):
            species_tag = ETree.SubElement(species_list,
                                           self._sbml_tag('species'))
            species_tag.set(self._sbml_tag('id'), 'M_' + species_id)
            species_tag.set(
                self._sbml_tag('name'),
                text_type(compound_name.get(species.name, species.name)))
            species_tag.set(
                self._sbml_tag('compartment'),
                model_compartments[species.compartment])
            species_tag.set(self._sbml_tag('constant'), 'false')
            species_tag.set(self._sbml_tag('boundaryCondition'), 'false')
            species_tag.set(self._sbml_tag('hasOnlySubstanceUnits'), 'true')
            if 'charge' in compound_properties[species.name]:
                species_tag.set(_tag('charge', FBC_V2), text_type(
                    compound_properties[species.name]['charge']))
            if 'formula' in compound_properties[species.name]:
                species_tag.set(_tag(
                    'chemicalFormula', FBC_V2), text_type(
                        compound_properties[species.name]['formula']))

            notes_tag = ETree.SubElement(species_tag, self._sbml_tag('notes'))
            body_tag = ETree.SubElement(notes_tag, _tag('body', XHTML_NS))
            self._add_properties_notes(
                body_tag, compound_properties[species.name])

        params_list = ETree.SubElement(
            model_tag, self._sbml_tag('listOfParameters'))

        # Create mapping for reactions containing flux limit definitions
        for rxn_id, lower_lim, upper_lim in itervalues(model.limits):
            flux_limits[rxn_id] = lower_lim, upper_lim

        params = {}
        gene_ids = {}

        if biomass_id is not None:
            self._add_fbc_objective(model_tag, biomass_id)

        # Create list of reactions
        reactions = ETree.SubElement(
            model_tag, self._sbml_tag('listOfReactions'))
        for eq_id, properties in sorted(iteritems(reaction_properties)):
            reaction_tag = ETree.SubElement(reactions,
                                            self._sbml_tag('reaction'))
            equation = properties['equation']

            reaction_tag.set(self._sbml_tag('id'), 'R_' + eq_id)
            if 'name' in properties:
                reaction_tag.set(self._sbml_tag('name'), properties['name'])
            reaction_tag.set(self._sbml_tag('reversible'), text_type(
                equation.direction == Direction.Both).lower())
            reaction_tag.set(self._sbml_tag('fast'), 'false')
            lower_str, upper_str = self._get_flux_bounds(
                eq_id, model, flux_limits, equation)

            params[upper_str] = 'P_'+self._make_safe_numerical_id(upper_str)
            params[lower_str] = 'P_'+self._make_safe_numerical_id(lower_str)
            reaction_tag.set(
                _tag('upperFluxBound', FBC_V2), params[upper_str])
            reaction_tag.set(
                _tag('lowerFluxBound', FBC_V2), params[lower_str])

            if 'genes' in properties:
                self._add_gene_associations(
                    eq_id, properties['genes'], gene_ids, reaction_tag)

            if any(value < 0 for _, value in equation.compounds):
                reactants = ETree.SubElement(
                    reaction_tag, self._sbml_tag('listOfReactants'))

            if any(value > 0 for _, value in equation.compounds):
                products = ETree.SubElement(
                    reaction_tag, self._sbml_tag('listOfProducts'))

            for compound, value in sorted(equation.compounds):
                dest_list = reactants if value < 0 else products
                spec_ref = ETree.SubElement(
                    dest_list, self._sbml_tag('speciesReference'))
                spec_ref.set(
                    self._sbml_tag('species'), 'M_' + model_species[compound])
                spec_ref.set(
                    self._sbml_tag('constant'), 'true')
                spec_ref.set(
                    self._sbml_tag('stoichiometry'), text_type(abs(value)))

            notes_tag = ETree.SubElement(reaction_tag, self._sbml_tag('notes'))
            body_tag = ETree.SubElement(notes_tag, _tag('body', XHTML_NS))
            self._add_properties_notes(body_tag, reaction_properties[eq_id])

            if self._cobra_flux_bounds is True:
                # Create COBRA-compliant parameter list
                kl_tag = ETree.SubElement(
                    reaction_tag, self._sbml_tag('kineticLaw'))
                math_tag = ETree.SubElement(kl_tag, self._sbml_tag('math'))
                ci_tag = ETree.SubElement(math_tag, _tag('ci', MATHML_NS))
                ci_tag.text = 'FLUX_VALUE'
                param_list = ETree.SubElement(
                    kl_tag, self._sbml_tag('listOfParameters'))

                ETree.SubElement(param_list, self._sbml_tag('parameter'), {
                    self._sbml_tag('id'): 'LOWER_BOUND',
                    self._sbml_tag('name'): 'LOWER_BOUND',
                    self._sbml_tag('value'): lower_str,
                    self._sbml_tag('constant'): 'true'
                })
                ETree.SubElement(param_list, self._sbml_tag('parameter'), {
                    self._sbml_tag('id'): 'UPPER_BOUND',
                    self._sbml_tag('name'): 'UPPER_BOUND',
                    self._sbml_tag('value'): upper_str,
                    self._sbml_tag('constant'): 'true'
                })

        for val, id in iteritems(params):
            param_tag = ETree.SubElement(
                params_list, self._sbml_tag('parameter'))
            param_tag.set(self._sbml_tag('id'), id)
            param_tag.set(self._sbml_tag('value'), val)
            param_tag.set(self._sbml_tag('constant'), 'true')

        self._add_gene_list(model_tag, gene_ids)

        tree = ETree.ElementTree(root)
        if pretty:
            self._indent(root)

        write_options = dict(
            encoding='ascii',
            default_namespace=self._namespace,
            xml_declaration=True
        )
        if PY3:
            write_options['encoding'] = 'unicode'

        tree.write(file, **write_options)


def convert_sbml_model(model):
    """Convert raw SBML model to extended model.

    Args:
        model: :class:`NativeModel` obtained from :class:`SBMLReader`.
    """
    biomass_reactions = set()
    for reaction in model.reactions:
        # Extract limits
        if reaction.id not in model.limits:
            lower, upper = parse_flux_bounds(reaction)
            if lower is not None or upper is not None:
                model.limits[reaction.id] = reaction.id, lower, upper

        # Detect objective
        objective = parse_objective_coefficient(reaction)
        if objective is not None and objective != 0:
            biomass_reactions.add(reaction.id)

    if len(biomass_reactions) == 1:
        model.biomass_reaction = next(iter(biomass_reactions))

    # Convert model to mutable entries
    convert_model_entries(model)

    # Detect extracelluar compartment
    if model.extracellular_compartment is None:
        extracellular = detect_extracellular_compartment(model)
        model.extracellular_compartment = extracellular

    # Convert exchange reactions to exchange compounds
    convert_exchange_to_compounds(model)


def entry_id_from_cobra_encoding(cobra_id):
    """Convert COBRA-encoded ID string to decoded ID string."""
    for escape, symbol in iteritems(_COBRA_DECODE_ESCAPES):
        cobra_id = cobra_id.replace(escape, symbol)
    return cobra_id


def create_convert_sbml_id_function(
        compartment_prefix='C_', reaction_prefix='R_',
        compound_prefix='M_', decode_id=entry_id_from_cobra_encoding):
    """Create function for converting SBML IDs.

    The returned function will strip prefixes, decode the ID using the provided
    function. These prefixes are common on IDs in SBML models because the IDs
    live in a global namespace.
    """
    def convert_sbml_id(entry):
        if isinstance(entry, BaseCompartmentEntry):
            prefix = compartment_prefix
        elif isinstance(entry, BaseReactionEntry):
            prefix = reaction_prefix
        elif isinstance(entry, BaseCompoundEntry):
            prefix = compound_prefix

        new_id = entry.id
        if decode_id is not None:
            new_id = decode_id(new_id)
        if prefix is not None and new_id.startswith(prefix):
            new_id = new_id[len(prefix):]

        return new_id

    return convert_sbml_id


def translate_sbml_compartment(entry, new_id):
    """Translate SBML compartment entry."""
    return DictCompartmentEntry(entry, id=new_id)


def translate_sbml_reaction(entry, new_id, compartment_map, compound_map):
    """Translate SBML reaction entry."""
    new_entry = DictReactionEntry(entry, id=new_id)

    # Convert compound IDs in reaction equation
    if new_entry.equation is not None:
        compounds = []
        for compound, value in new_entry.equation.compounds:
            # Translate compartment to new ID, if available.
            compartment = compartment_map.get(
                compound.compartment, compound.compartment)
            new_compound = compound.translate(
                lambda name: compound_map.get(name, name)).in_compartment(
                    compartment)
            compounds.append((new_compound, value))

        new_entry.equation = Reaction(
            new_entry.equation.direction, compounds)

    # Get XHTML notes properties
    for key, value in iteritems(parse_xhtml_reaction_notes(entry)):
        if key not in new_entry.properties:
            new_entry.properties[key] = value

    return new_entry


def translate_sbml_compound(entry, new_id, compartment_map):
    """Translate SBML compound entry."""
    new_entry = DictCompoundEntry(entry, id=new_id)

    if 'compartment' in new_entry.properties:
        old_compartment = new_entry.properties['compartment']
        new_entry.properties['compartment'] = compartment_map.get(
            old_compartment, old_compartment)

    # Get XHTML notes properties
    for key, value in iteritems(parse_xhtml_species_notes(entry)):
        if key not in new_entry.properties:
            new_entry.properties[key] = value

    return new_entry


def convert_model_entries(
        model, convert_id=create_convert_sbml_id_function(),
        create_unique_id=None,
        translate_compartment=translate_sbml_compartment,
        translate_reaction=translate_sbml_reaction,
        translate_compound=translate_sbml_compound):
    """Convert and decode model entries.

    Model entries are converted to new entries using the translate functions
    and IDs are converted using the given coversion function. If ID conversion
    would create a clash of IDs, the ``create_unique_id`` function is called
    with a container of current IDs and the base ID to generate a unique ID
    from. The translation functions take an existing entry and the new ID.

    All references within the model are updated to use new IDs: compartment
    boundaries, limits, exchange, model, biomass reaction, etc.

    Args:
        model: :class:`NativeModel`.
    """
    def find_new_ids(entries):
        """Create new IDs for entries."""
        id_map = {}
        new_ids = set()
        for entry in entries:
            new_id = convert_id(entry)
            if new_id in new_ids:
                if create_unique_id is not None:
                    new_id = create_unique_id(new_ids, new_id)
                else:
                    raise ValueError(
                        'Entity ID {!r} is not unique after conversion'.format(
                            entry.id))

            id_map[entry.id] = new_id
            new_ids.add(new_id)

        return id_map

    # Find new IDs for all entries
    compartment_map = find_new_ids(model.compartments)
    compound_map = find_new_ids(model.compounds)
    reaction_map = find_new_ids(model.reactions)

    # Create new compartment entries
    new_compartments = []
    for compartment in model.compartments:
        new_id = compartment_map[compartment.id]
        new_compartments.append(translate_compartment(compartment, new_id))

    # Create new compound entries
    new_compounds = []
    for compound in model.compounds:
        new_id = compound_map[compound.id]
        new_compounds.append(
            translate_compound(compound, new_id, compartment_map))

    # Create new reaction entries
    new_reactions = []
    for reaction in model.reactions:
        new_id = reaction_map[reaction.id]
        new_entry = translate_reaction(
            reaction, new_id, compartment_map, compound_map)
        new_reactions.append(new_entry)

    # Update entries
    model.compartments.clear()
    model.compartments.update(new_compartments)

    model.compounds.clear()
    model.compounds.update(new_compounds)

    model.reactions.clear()
    model.reactions.update(new_reactions)

    # Convert compartment boundaries
    new_boundaries = []
    for boundary in model.compartment_boundaries:
        c1, c2 = (compartment_map.get(c, c) for c in boundary)
        new_boundaries.append(tuple(sorted(c1, c2)))

    model.compartment_boundaries.clear()
    model.compartment_boundaries.update(new_boundaries)

    # Convert limits
    new_limits = []
    for reaction, lower, upper in itervalues(model.limits):
        new_reaction_id = reaction_map.get(reaction, reaction)
        new_limits.append((new_reaction_id, lower, upper))

    model.limits.clear()
    model.limits.update((limit[0], limit) for limit in new_limits)

    # Convert exchange
    new_exchanges = []
    for compound, reaction, lower, upper in itervalues(model.exchange):
        new_compound_id = compound.translated(
            lambda name: compound_map.get(name, name))
        new_reaction_id = reaction_map.get(reaction, reaction)
        new_exchanges.append((new_compound_id, new_reaction_id, lower, upper))

    model.exchange.clear()
    model.exchange.update((ex[0], ex) for ex in new_exchanges)

    # Convert model
    new_model = []
    for reaction in model.model:
        new_id = reaction_map.get(reaction, reaction)
        new_model.append(new_id)

    model.model.clear()
    model.model.update((new_id, None) for new_id in new_model)

    # Convert other properties
    if model.biomass_reaction is not None:
        old_id = model.biomass_reaction
        model.biomass_reaction = reaction_map.get(old_id, old_id)

    if model.extracellular_compartment is not None:
        old_id = model.extracellular_compartment
        model.extracellular_compartment = compartment_map.get(old_id, old_id)

    if model.default_compartment is not None:
        old_id = model.default_compartment
        model.default_compartment = compartment_map.get(old_id, old_id)


def parse_xhtml_notes(entry):
    """Yield key, value pairs parsed from the XHTML notes section.

    Each key, value pair must be defined in its own text block, e.g.
    ``<p>key: value</p><p>key2: value2</p>``. The key and value must be
    separated by a colon. Whitespace is stripped from both key and value, and
    quotes are removed from values if present. The key is normalized by
    conversion to lower case and spaces replaced with underscores.

    Args:
        entry: :class:`_SBMLEntry`.
    """
    for note in entry.xml_notes.itertext():
        m = re.match(r'^([^:]+):(.+)$', note)
        if m:
            key, value = m.groups()
            key = key.strip().lower().replace(' ', '_')
            value = value.strip()
            m = re.match(r'^"(.*)"$', value)
            if m:
                value = m.group(1)
            if value != '':
                yield key, value


def parse_xhtml_species_notes(entry):
    """Return species properties defined in the XHTML notes.

    Older SBML models often define additional properties in the XHTML notes
    section because structured methods for defining properties had not been
    developed. This will try to parse the following properties: ``PUBCHEM ID``,
    ``CHEBI ID``, ``FORMULA``, ``KEGG ID``, ``CHARGE``.

    Args:
        entry: :class:`SBMLSpeciesEntry`.
    """
    properties = {}
    if entry.xml_notes is not None:
        cobra_notes = dict(parse_xhtml_notes(entry))

        for key in ('pubchem_id', 'chebi_id'):
            if key in cobra_notes:
                properties[key] = cobra_notes[key]

        if 'formula' in cobra_notes:
            properties['formula'] = cobra_notes['formula']

        if 'kegg_id' in cobra_notes:
            properties['kegg'] = cobra_notes['kegg_id']

        if 'charge' in cobra_notes:
            try:
                value = int(cobra_notes['charge'])
            except ValueError:
                logger.warning(
                    'Unable to parse charge for {} as an'
                    ' integer: {}'.format(
                        entry.id, cobra_notes['charge']))
                value = cobra_notes['charge']
            properties['charge'] = value

    return properties


def parse_xhtml_reaction_notes(entry):
    """Return reaction properties defined in the XHTML notes.

    Older SBML models often define additional properties in the XHTML notes
    section because structured methods for defining properties had not been
    developed. This will try to parse the following properties: ``SUBSYSTEM``,
    ``GENE ASSOCIATION``, ``EC NUMBER``, ``AUTHORS``, ``CONFIDENCE``.

    Args:
        entry: :class:`SBMLReactionEntry`.
    """
    properties = {}
    if entry.xml_notes is not None:
        cobra_notes = dict(parse_xhtml_notes(entry))

        if 'subsystem' in cobra_notes:
            properties['subsystem'] = cobra_notes['subsystem']

        if 'gene_association' in cobra_notes:
            properties['genes'] = cobra_notes['gene_association']

        if 'ec_number' in cobra_notes:
            properties['ec'] = cobra_notes['ec_number']

        if 'authors' in cobra_notes:
            properties['authors'] = [
                a.strip() for a in cobra_notes['authors'].split(';')]

        if 'confidence' in cobra_notes:
            try:
                value = int(cobra_notes['confidence'])
            except ValueError:
                logger.warning(
                    'Unable to parse confidence level for {} as an'
                    ' integer: {}'.format(
                        entry.id, cobra_notes['confidence']))
                value = cobra_notes['confidence']
            properties['confidence'] = value

    return properties


def parse_objective_coefficient(entry):
    """Return objective value for reaction entry.

    Detect objectives that are specified using the non-standardized
    kinetic law parameters which are used by many pre-FBC SBML models. The
    objective coefficient is returned for the given reaction, or None if
    undefined.

    Args:
        entry: :class:`SBMLReactionEntry`.
    """
    for parameter in entry.kinetic_law_reaction_parameters:
        pid, name, value, units = parameter
        if (pid == 'OBJECTIVE_COEFFICIENT' or
                name == 'OBJECTIVE_COEFFICIENT'):
            return value

    return None


def parse_flux_bounds(entry):
    """Return flux bounds for reaction entry.

    Detect flux bounds that are specified using the non-standardized
    kinetic law parameters which are used by many pre-FBC SBML models. The
    flux bounds are returned as a pair of lower, upper bounds. The returned
    bound is None if undefined.

    Args:
        entry: :class:`SBMLReactionEntry`.
    """
    lower_bound = None
    upper_bound = None
    for parameter in entry.kinetic_law_reaction_parameters:
        pid, name, value, units = parameter
        if pid == 'UPPER_BOUND' or name == 'UPPER_BOUND':
            upper_bound = value
        elif pid == 'LOWER_BOUND' or name == 'LOWER_BOUND':
            lower_bound = value

    return lower_bound, upper_bound


def detect_extracellular_compartment(model):
    """Detect the identifier for equations with extracellular compartments.

    Args:
        model: :class:`NativeModel`.
    """
    extracellular_key = Counter()

    for reaction in model.reactions:
        equation = reaction.equation
        if equation is None:
            continue

        if len(equation.compounds) == 1:
            compound, _ = equation.compounds[0]
            compartment = compound.compartment
            extracellular_key[compartment] += 1
    if len(extracellular_key) == 0:
        return None
    else:
        best_key, _ = extracellular_key.most_common(1)[0]

    logger.info('{} is extracellular compartment'.format(best_key))

    return best_key


def convert_exchange_to_compounds(model):
    """Convert exchange reactions in model to exchange compounds.

    Only exchange reactions in the extracellular compartment are converted.
    The extracelluar compartment must be defined for the model.

    Args:
        model: :class:`NativeModel`.
    """
    # Build set of exchange reactions
    exchanges = set()
    for reaction in model.reactions:
        equation = reaction.properties.get('equation')
        if equation is None:
            continue

        if len(equation.compounds) != 1:
            # Provide warning for exchange reactions with more than
            # one compound, they won't be put into the exchange definition
            if (len(equation.left) == 0) != (len(equation.right) == 0):
                logger.warning('Exchange reaction {} has more than one'
                               ' compound, it was not converted to'
                               ' exchange compound'.format(reaction.id))
            continue

        exchanges.add(reaction.id)

    # Convert exchange reactions into exchange compounds
    for reaction_id in exchanges:
        equation = model.reactions[reaction_id].equation
        compound, value = equation.compounds[0]
        if compound.compartment != model.extracellular_compartment:
            continue

        if compound in model.exchange:
            logger.warning(
                'Compound {} is already defined in the exchange'
                ' definition'.format(compound))
            continue

        # We multiply the flux bounds by value in order to create equivalent
        # exchange reactions with stoichiometric value of one. If the flux
        # bounds are not set but the reaction is unidirectional, the implicit
        # flux bounds must be used.
        lower_flux, upper_flux = None, None
        if reaction_id in model.limits:
            _, lower, upper = model.limits[reaction_id]
            if lower is not None:
                lower_flux = lower * abs(value)
            if upper is not None:
                upper_flux = upper * abs(value)

        if lower_flux is None and equation.direction == Direction.Forward:
            lower_flux = 0
        if upper_flux is None and equation.direction == Direction.Reverse:
            upper_flux = 0

        # If the stoichiometric value of the reaction is reversed, the flux
        # limits must be flipped.
        if value > 0:
            lower_flux, upper_flux = (
                -upper_flux if upper_flux is not None else None,
                -lower_flux if lower_flux is not None else None)

        model.exchange[compound] = (
            compound, reaction_id, lower_flux, upper_flux)

        model.reactions.discard(reaction_id)
        model.limits.pop(reaction_id, None)


def merge_equivalent_compounds(model):
    """Merge equivalent compounds in various compartments.

    Tries to detect and merge compound entries that represent the same
    compound in different compartments. The entries are only merged if all
    properties are equivalent. Compound entries must have an ID with a suffix
    of an underscore followed by the compartment ID. This suffix will be
    stripped and compounds with identical IDs are merged if the properties
    are identical.

    Args:
        model: :class:`NativeModel`.
    """
    def dicts_are_compatible(d1, d2):
        return all(key not in d1 or key not in d2 or d1[key] == d2[key]
                   for key in set(d1) | set(d2))

    compound_compartment = {}
    inelegible = set()
    for reaction in model.reactions:
        equation = reaction.equation
        if equation is None:
            continue

        for compound, _ in equation.compounds:
            compartment = compound.compartment
            if compartment is not None:
                compound_compartment[compound.name] = compartment
                if not compound.name.endswith('_{}'.format(compartment)):
                    inelegible.add(compound.name)

    compound_groups = {}
    for compound_id, compartment in iteritems(compound_compartment):
        if compound_id in inelegible:
            continue

        suffix = '_{}'.format(compound_compartment[compound_id])
        if compound_id.endswith(suffix):
            group_name = compound_id[:-len(suffix)]
            compound_groups.setdefault(group_name, set()).add(compound_id)

    compound_mapping = {}
    merged_compounds = {}
    for group, compound_set in iteritems(compound_groups):
        # Try to merge as many compounds as possible
        merged = []
        for compound_id in compound_set:
            props = dict(model.compounds[compound_id].properties)

            # Ignore differences in ID and compartment properties
            props.pop('id', None)
            props.pop('compartment', None)

            for merged_props, merged_set in merged:
                if dicts_are_compatible(props, merged_props):
                    merged_set.add(compound_id)
                    merged_props.update(props)
                    break
                else:
                    keys = set(key for key in set(props) | set(merged_props)
                               if key not in props or
                               key not in merged_props or
                               props[key] != merged_props[key])
                    logger.info(
                        'Unable to merge {} into {}, difference in'
                        ' keys: {}'.format(
                            compound_id, ', '.join(merged_set),
                            ', '.join(keys)))
            else:
                merged.append((props, {compound_id}))

        if len(merged) == 1:
            # Merge into one set with the group name
            merged_props, merged_set = merged[0]

            for compound_id in merged_set:
                compound_mapping[compound_id] = group
            merged_compounds[group] = merged_props
        else:
            # Since we cannot merge all compounds, create new group names
            # based on the group and compartments.
            for merged_props, merged_set in merged:
                compartments = set(compound_compartment[c] for c in merged_set)
                merged_name = '{}_{}'.format(
                    group, '_'.join(sorted(compartments)))

                for compound_id in merged_set:
                    compound_mapping[compound_id] = merged_name
                merged_compounds[merged_name] = merged_props

    # Translate reaction compounds
    for reaction in model.reactions:
        equation = reaction.equation
        if equation is None:
            continue

        reaction.equation = equation.translated_compounds(
            lambda c: compound_mapping.get(c, c))

    # Translate compound entries
    new_compounds = []
    for compound in model.compounds:
        if compound.id not in compound_mapping:
            new_compounds.append(compound)
        else:
            group = compound_mapping[compound.id]
            if group not in merged_compounds:
                continue
            props = merged_compounds.pop(group)
            props['id'] = group
            new_compounds.append(DictCompoundEntry(
                props, filemark=compound.filemark))

    model.compounds.clear()
    model.compounds.update(new_compounds)

    # Translate exchange
    new_exchange = OrderedDict()
    for compound, reaction_id, lower, upper in itervalues(model.exchange):
        new_compound = compound.translate(
            lambda name: compound_mapping.get(name, name))
        new_exchange[new_compound] = new_compound, reaction_id, lower, upper

    model.exchange.clear()
    model.exchange.update(new_exchange)
