
"""Parser for SBML model files"""

import xml.etree.ElementTree as ET
from decimal import Decimal
from fractions import Fraction
from functools import partial
from itertools import count

from ..database import MetabolicDatabase, DictDatabase
from ..reaction import Reaction, Compound


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


def _tag(tag, namespace=None):
    """Prepend namespace to tag name"""
    if namespace is None:
        return str(tag)
    return '{{{}}}{}'.format(namespace, tag)


class ParseError(Exception):
    """Error parsing SBML file"""


class SBMLDatabase(MetabolicDatabase):
    """Reader of SBML model files

    The constructor takes a file-like object which will be parsed as XML and
    then as SBML according to the specification. If the ``strict`` parameter is
    set to False, the parser will revert to a more lenient parsing which is
    required for many older models. This tries to mimic the inconsistencies
    employed by COBRA when parsing models.
    """

    def __init__(self, file, strict=False):
        # Parse SBML file
        tree = ET.parse(file)
        root = tree.getroot()

        self._strict = strict

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
        self._database = DictDatabase()

        # Compounds
        self._model_compounds = {}
        self._species = self._model.find(self._sbml_tag('listOfSpecies'))
        for species in self._species.iterfind(self._sbml_tag('species')):
            species_name = species.get('name')
            species_id = self._element_get_id(species)
            species_comp = species.get('compartment')
            self._model_compounds[species_id] = species_name, species_comp

            species_boundary = species.get('boundaryCondition', False)

            # In non-strict mode the species that ends with _b are considered
            # boundary conditions.
            if not self._strict and species_id.endswith('_b'):
                species_boundary = True

            # Add implicit exchange reaction if compound is boundary condition
            if species_boundary:
                reaction_name = species_id + '_impl_EX'
                reaction = Reaction(
                    Reaction.Bidir,
                    [(Compound(species_id, compartment=species_comp), 1)], [])
                self._database.set_reaction(reaction_name, reaction)

        # Reactions
        self._reactions = self._model.find(self._sbml_tag('listOfReactions'))
        for reaction in self._reactions.iterfind(self._sbml_tag('reaction')):
            reaction_name = reaction.get('name')
            reaction_id = self._element_get_id(reaction)
            reaction_rev = reaction.get('reversible', 'true') == 'true'

            left = []
            for species_id, value in self._parse_species_references(
                    reaction, 'listOfReactants'):
                if species_id not in self._model_compounds:
                    if not self._strict:
                        # In non-strict mode simply skip these references
                        continue
                    raise ParseError(
                        'Reaction {} references non-existent'
                        ' species {}'.format(reaction_id, species_id))

                species_name, species_comp = self._model_compounds[species_id]
                compound = Compound(species_id, compartment=species_comp)
                left.append((compound, value))

            right = []
            for species_id, value in self._parse_species_references(
                    reaction, 'listOfProducts'):
                if species_id not in self._model_compounds:
                    if not self._strict:
                        # In non-strict mode simply skip these references
                        continue
                    raise ParseError(
                        'Reaction {} references non-existent'
                        ' species {}'.format(reaction_id, species_id))

                species_name, species_comp = self._model_compounds[species_id]
                compound = Compound(species_id, compartment=species_comp)
                right.append((compound, value))

            # Add reaction to database
            direction = Reaction.Bidir if reaction_rev else Reaction.Right
            self._database.set_reaction(reaction_id, Reaction(direction, left, right))

    def _element_get_id(self, element):
        """Get id of reaction or species element

        In old levels the name is used as the id. This method returns the correct
        attribute depending on the level.
        """

        return element.get('id') if self._level > 1 else element.get('name')

    def _parse_species_references(self, root, name):
        """Yield species id and parsed value for a speciesReference list"""
        for species in root.iterfind('./{}/{}'.format(
                self._sbml_tag(name), self._sbml_tag('speciesReference'))):

            species_id = species.get('species')

            if self._level == 1:
                # In SBML level 1 only positive integers are allowed for
                # stoichiometry but a positive integer denominator can be
                # specified.
                try:
                    value = int(species.get('stoichiometry', 1))
                    denom = int(species.get('denominator', 1))
                    species_value = Fraction(value, denom)
                except ValueError:
                    if self._strict:
                        raise
                    species_value = Decimal(species.get('stoichiometry', 1))
            elif self._level == 2:
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
            elif self._level == 3:
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
    def reactions(self):
        return self._database.reactions

    @property
    def compounds(self):
        return self._database.compounds

    def has_reaction(self, reaction_id):
        return self._database.has_reaction(reaction_id)

    def is_reversible(self, reaction_id):
        return self._database.is_reversible(reaction_id)

    def get_reaction_values(self, reaction_id):
        return self._database.get_reaction_values(reaction_id)

    def get_compound_reactions(self, compound):
        return self._database.get_compound_reactions(compound)

    def get_compound_name(self, compound):
        """Name of compound"""
        if compound.name not in self._model_compounds:
            raise ValueError('Unknown compound: {}'.format(compound))
        name, comp = self._model_compounds[compound.name]
        return name

    def get_reaction_notes_elements(self):
        """Yield tuples of reaction ids, and notes as ElementTree elements"""
        for reaction in self._reactions.iterfind(self._sbml_tag('reaction')):
            reaction_id = self._element_get_id(reaction)
            yield reaction_id, reaction.find(self._sbml_tag('notes'))

    def get_compound_notes_elements(self):
        """Yield tuples of compound ids, and notes as ElementTree elements"""
        for compound in self._species.iterfind(self._sbml_tag('species')):
            compound_id = self._element_get_id(compound)
            yield compound_id, compound.find(self._sbml_tag('notes'))

    def get_kinetic_law_reaction_parameters(self):
        """Yield tuples of reaction ids and the value of a kinetic law reaction
        parameters
        """
        for reaction in self._reactions.iterfind(self._sbml_tag('reaction')):
            reaction_id = self._element_get_id(reaction)
            for parameter in reaction.iterfind(
                    './{}/{}/{}'.format(self._sbml_tag('kineticLaw'),
                                        self._sbml_tag('listOfParameters'),
                                        self._sbml_tag('parameter'))):
                param_id = parameter.get('id')
                param_value = float(parameter.get('value'))
                param_units = parameter.get('units')

                yield reaction_id, (param_id, param_value, param_units)


class SBMLWriter(object):
    """Writer of SBML files"""

    def __init__(self):
        self._namespace = SBML_NS_L3_V1_CORE
        self._sbml_tag = partial(_tag, namespace=self._namespace)

    def write_model(self, file, model, compounds):
        """Write a given model to file"""

        # Load compound information
        compound_name = {}
        for compound in compounds:
            compound_name[compound.id] = compound.name if compound.name is not None else compound.id

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
                    model_compartments[compound.compartment] = next(compartment_id)
                if compound not in model_species:
                    model_species[compound] = next(compound_id)

        # Create list of compartments
        compartments = ET.SubElement(model_tag, self._sbml_tag('listOfCompartments'))
        for compartment, compartment_id in model_compartments.iteritems():
            compartment_tag = ET.SubElement(compartments, self._sbml_tag('compartment'))
            compartment_tag.set(self._sbml_tag('id'), compartment_id)
            compartment_tag.set(self._sbml_tag('name'), str(compartment))

        # Create list of species
        species_list = ET.SubElement(model_tag, self._sbml_tag('listOfSpecies'))
        for species, species_id in model_species.iteritems():
            species_tag = ET.SubElement(species_list, self._sbml_tag('species'))
            species_tag.set(self._sbml_tag('id'), species_id)
            species_tag.set(self._sbml_tag('name'), str(species.translate(lambda x: compound_name.get(x, x))))
            species_tag.set(self._sbml_tag('compartment'), model_compartments[species.compartment])

        # Create list of reactions
        reactions = ET.SubElement(model_tag, self._sbml_tag('listOfReactions'))
        for reaction in model.reactions:
            reaction_tag = ET.SubElement(reactions, self._sbml_tag('reaction'))
            reaction_tag.set(self._sbml_tag('id'), next(reaction_id))
            reaction_tag.set(self._sbml_tag('name'), reaction)
            reaction_tag.set(self._sbml_tag('reversible'), 'true' if model.is_reversible(reaction) else 'false')

            reactants = ET.SubElement(reaction_tag, self._sbml_tag('listOfReactants'))
            products = ET.SubElement(reaction_tag, self._sbml_tag('listOfProducts'))

            for compound, value in model.get_reaction_values(reaction):
                dest_list = reactants if value < 0 else products
                spec_ref = ET.SubElement(dest_list, self._sbml_tag('speciesReference'))
                spec_ref.set(self._sbml_tag('species'), model_species[compound])
                spec_ref.set(self._sbml_tag('stoichiometry'), str(value))

        tree = ET.ElementTree(root)
        tree.write(file, default_namespace=self._namespace)
