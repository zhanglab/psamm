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

"""Parser for SBML model files"""

from __future__ import unicode_literals

import xml.etree.ElementTree as ET
from decimal import Decimal
from fractions import Fraction
from functools import partial
from itertools import count
import logging

from six import itervalues, iteritems

from ..reaction import Reaction, Compound, Direction


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


def _tag(tag, namespace=None):
    """Prepend namespace to tag name"""
    if namespace is None:
        return str(tag)
    return '{{{}}}{}'.format(namespace, tag)


class ParseError(Exception):
    """Error parsing SBML file"""


class _SBMLEntry(object):
    """Base class for compound and reaction entries"""

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
        """Access the entity notes as an XML document fragment"""
        return self._root.find(self._reader._sbml_tag('notes'))


class SpeciesEntry(_SBMLEntry):
    """Species entry in the SBML file"""

    def __init__(self, reader, root):
        super(SpeciesEntry, self).__init__(reader, root)

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
        if 'charge' in self._root.attrib and self._root.get('charge') != '':
            properties['charge'] = int(self._root.get('charge'))

        return properties


class ReactionEntry(_SBMLEntry):
    """Reaction entry in SBML file"""

    def __init__(self, reader, root):
        super(ReactionEntry, self).__init__(reader, root)

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
            param_value = float(parameter.get('value'))
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


class ObjectiveEntry(object):
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
            raise ParseError('Missing "type" attribute on objective: {}'.formt(
                self.id))

        # Find flux objectives
        self._reactions = {}
        for fo in self._root.iterfind('./{}/{}'.format(
                _tag('listOfFluxObjectives', self._namespace),
                _tag('fluxObjective', self._namespace))):
            reaction = fo.get(_tag('reaction', self._namespace))
            coefficient = float(fo.get(_tag('coefficient', self._namespace)))
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


class FluxBoundEntry(object):
    """Flux bound defined with FBC"""

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
        elif self._operation not in ('lessEqual', 'greaterEqual', 'equal'):
            raise ParseError('Invalid flux bound operation: {}'.format(
                self._operation))

        value = self._root.get(_tag('value', self._namespace))
        if value is None:
            raise ParseError('Flux bound is missing "value" attribute')
        try:
            self._value = float(value)
        except ValueError:
            raise ParseError('Unable to parse flux bound value')

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        return self._name

    @property
    def reaction(self):
        return self._reaction

    @property
    def operation(self):
        return self._operation

    @property
    def value(self):
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
    the reaction equations.
    """

    def __init__(self, file, strict=False, ignore_boundary=False):
        # Parse SBML file
        tree = ET.parse(file)
        root = tree.getroot()

        self._strict = strict
        self._ignore_boundary = ignore_boundary

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
                    value = float(param.get('value'))
                    self._model_constants[param_id] = value

        # Flux bounds
        self._flux_bounds = []
        self._reaction_flux_bounds = {}
        if self._level == 3:
            # Only suported in FBC V1
            flux_bounds = self._model.find(_tag('listOfFluxBounds', FBC_V1))
            if flux_bounds is not None:
                for flux_bound in flux_bounds.iterfind(
                        _tag('fluxBound', FBC_V1)):
                    entry = FluxBoundEntry(self, FBC_V1, flux_bound)
                    self._flux_bounds.append(entry)

                    # Create reference from reaction to flux bound
                    entries = self._reaction_flux_bounds.setdefault(
                        entry.reaction, [])
                    entries.append(entry)

        # Species
        self._model_species = {}
        self._species = self._model.find(self._sbml_tag('listOfSpecies'))
        for species in self._species.iterfind(self._sbml_tag('species')):
            entry = SpeciesEntry(self, species)
            self._model_species[entry.id] = entry

        # Reactions
        self._model_reactions = {}
        self._reactions = self._model.find(self._sbml_tag('listOfReactions'))
        for reaction in self._reactions.iterfind(self._sbml_tag('reaction')):
            entry = ReactionEntry(self, reaction)
            self._model_reactions[entry.id] = entry

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
                    entry = ObjectiveEntry(self, fbc_ns, objective)
                    self._model_objectives[entry.id] = entry

                active = objectives.get(_tag('activeObjective', fbc_ns))
                if active is None or active not in self._model_objectives:
                    raise ParseError('Active objective is invalid: {}'.format(
                        active))

                self._active_objective = self._model_objectives[active]

    def get_reaction(self, reaction_id):
        """Return :class:`.ReactionEntry` corresponding to reaction_id"""
        return self._model_reactions[reaction_id]

    def get_species(self, species_id):
        """Return :class:`.SpeciesEntry` corresponding to species_id"""
        return self._model_species[species_id]

    def get_objective(self, objective_id):
        """Return :class:`.ObjectiveEntry` corresponding to objective_id"""
        return self._model_objectives[objective_id]

    def get_active_objective(self):
        return self._active_objective

    @property
    def reactions(self):
        """Iterator over :class:`ReactionEntries <.ReactionEntry>`"""
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
        """Iterator over :class:`.ObjectiveEntry`"""
        return itervalues(self._model_objectives)

    @property
    def flux_bounds(self):
        """Iterator over :class:`.FluxBoundEntry`"""
        return iter(self._flux_bounds)

    @property
    def id(self):
        """Model ID"""
        return self._model.get('id', None)

    @property
    def name(self):
        """Model name"""
        return self._model.get('name', None)


class SBMLWriter(object):
    """Writer of SBML files"""

    def __init__(self):
        self._namespace = SBML_NS_L3_V1_CORE
        self._sbml_tag = partial(_tag, namespace=self._namespace)

    def write_model(self, file, model, reactions, compounds):
        """Write a given model to file"""

        ET.register_namespace('mathml', MATHML_NS)
        ET.register_namespace('xhtml', XHTML_NS)

        # Load compound information
        compound_name = {}
        for compound in compounds:
            compound_name[compound.id] = (
                compound.name if compound.name is not None else compound.id)

        reaction_genes = {}
        for reaction in reactions:
            if reaction.genes is not None:
                reaction_genes[reaction.id] = reaction.genes

        root = ET.Element(self._sbml_tag('sbml'))
        root.set(self._sbml_tag('level'), '3')
        root.set(self._sbml_tag('version'), '1')

        model_tag = ET.SubElement(root, self._sbml_tag('model'))

        # Generators of unique IDs
        compound_id = ('M_'+str(i) for i in count(1))
        compartment_id = ('C_'+str(i) for i in count(1))
        reaction_id = ('R_'+str(i) for i in count(1))

        # Build mapping from Compound to species ID
        model_compartments = {}
        model_species = {}
        for reaction in model.reactions:
            for compound, value in model.get_reaction_values(reaction):
                if compound.compartment not in model_compartments:
                    model_compartments[compound.compartment] = (
                        next(compartment_id))
                if compound not in model_species:
                    model_species[compound] = next(compound_id)

        # Create list of compartments
        compartments = ET.SubElement(
            model_tag, self._sbml_tag('listOfCompartments'))
        for compartment, compartment_id in iteritems(model_compartments):
            compartment_tag = ET.SubElement(
                compartments, self._sbml_tag('compartment'))
            compartment_tag.set(self._sbml_tag('id'), compartment_id)
            compartment_tag.set(self._sbml_tag('name'), str(compartment))

        # Create list of species
        species_list = ET.SubElement(
            model_tag, self._sbml_tag('listOfSpecies'))
        for species, species_id in iteritems(model_species):
            species_tag = ET.SubElement(species_list,
                                        self._sbml_tag('species'))
            species_tag.set(self._sbml_tag('id'), species_id)
            species_tag.set(
                self._sbml_tag('name'),
                str(species.translate(lambda x: compound_name.get(x, x))))
            species_tag.set(
                self._sbml_tag('compartment'),
                model_compartments[species.compartment])

        # Create list of reactions
        reactions = ET.SubElement(model_tag, self._sbml_tag('listOfReactions'))
        for reaction in model.reactions:
            reaction_tag = ET.SubElement(reactions, self._sbml_tag('reaction'))
            reaction_tag.set(self._sbml_tag('id'), next(reaction_id))
            reaction_tag.set(self._sbml_tag('name'), reaction)
            reaction_tag.set(
                self._sbml_tag('reversible'),
                'true' if model.is_reversible(reaction) else 'false')

            reactants = ET.SubElement(
                reaction_tag, self._sbml_tag('listOfReactants'))
            products = ET.SubElement(
                reaction_tag, self._sbml_tag('listOfProducts'))

            for compound, value in model.get_reaction_values(reaction):
                dest_list = reactants if value < 0 else products
                spec_ref = ET.SubElement(
                    dest_list, self._sbml_tag('speciesReference'))
                spec_ref.set(
                    self._sbml_tag('species'), model_species[compound])
                spec_ref.set(self._sbml_tag('stoichiometry'), str(abs(value)))

            notes_tag = ET.SubElement(reaction_tag, self._sbml_tag('notes'))
            if reaction in reaction_genes:
                body_tag = ET.SubElement(notes_tag, _tag('body', XHTML_NS))
                p_tag = ET.SubElement(body_tag, _tag('p', XHTML_NS))
                p_tag.text = 'GENE_ASSOCIATION: {}'.format(
                    reaction_genes[reaction])

            # Create COBRA-compliant parameter list
            kl_tag = ET.SubElement(reaction_tag, self._sbml_tag('kineticLaw'))
            math_tag = ET.SubElement(kl_tag, self._sbml_tag('math'))
            ci_tag = ET.SubElement(math_tag, _tag('ci', MATHML_NS))
            ci_tag.text = 'FLUX_VALUE'
            param_list = ET.SubElement(
                kl_tag, self._sbml_tag('listOfParameters'))
            ET.SubElement(param_list, self._sbml_tag('parameter'), {
                self._sbml_tag('id'): 'LOWER_BOUND',
                self._sbml_tag('name'): 'LOWER_BOUND',
                self._sbml_tag('value'): str(model.limits[reaction].lower)
            })
            ET.SubElement(param_list, self._sbml_tag('parameter'), {
                self._sbml_tag('id'): 'UPPER_BOUND',
                self._sbml_tag('name'): 'UPPER_BOUND',
                self._sbml_tag('value'): str(model.limits[reaction].upper)
            })

        tree = ET.ElementTree(root)
        tree.write(file, default_namespace=self._namespace)
