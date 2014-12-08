
'''Parser for SBML model files'''

import xml.etree.ElementTree as ET
from decimal import Decimal
from functools import partial

from .database import MetabolicDatabase, DictDatabase
from .reaction import Reaction, Compound

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

def tag(tag, namespace=None):
    '''Prepend namespace to tag name'''
    if namespace is None:
        return str(tag)
    return '{{{}}}{}'.format(namespace, tag)

class ParseError(Exception):
    '''Error parsing SBML file'''

class SBMLDatabase(MetabolicDatabase):
    '''Reaction database backed by an SBML file'''

    def __init__(self, file):
        # Parse SBML file
        tree = ET.parse(file)
        root = tree.getroot()

        # Parse level and version
        self._sbml_tag = None
        self._level = int(root.get('level'))
        self._version = int(root.get('version'))

        if self._level == 1:
            self._sbml_tag = partial(tag, namespace=SBML_NS_L1)
        elif self._level == 2:
            if self._version == 1:
                self._sbml_tag = partial(tag, namespace=SBML_NS_L2)
            elif self._version == 2:
                self._sbml_tag = partial(tag, namespace=SBML_NS_L2_V2)
            elif self._version == 3:
                self._sbml_tag = partial(tag, namespace=SBML_NS_L2_V3)
            elif self._version == 4:
                self._sbml_tag = partial(tag, namespace=SBML_NS_L2_V4)
            elif self._version == 5:
                self._sbml_tag = partial(tag, namespace=SBML_NS_L2_V5)
            else:
                raise ParseError('SBML level 2, version {} not implemented'.format(self._version))
        elif self._level == 3:
            if self._version == 1:
                self._sbml_tag = partial(tag, namespace=SBML_NS_L3_V1_CORE)
            else:
                raise ParseError('SBML level 3, version {} not implemented'.format(self._version))
        else:
            raise ParseError('SBML level {} not implemented'.format(self._level))

        self._model = root.find(self._sbml_tag('model'))
        self._database = DictDatabase()

        # Compounds
        self._model_compounds = {}
        self._species = self._model.find(self._sbml_tag('listOfSpecies'))
        for species in self._species.iterfind(self._sbml_tag('species')):
            species_name = species.get('name')
            species_id = species.get('id') if self._level > 1 else species_name
            species_comp = species.get('compartment')
            self._model_compounds[species_id] = species_name, species_comp

            # Add implicit exchange reaction if compound is boundary condition
            species_boundary = species.get('boundaryCondition', False)
            if species_boundary:
                reaction_name = species_id + '_impl_EX'
                reaction = Reaction(Reaction.Bidir, [(Compound(species_id, compartment=species_comp), 1)], [])
                self._database.set_reaction(reaction_name, reaction)

        # Reactions
        self._reactions = self._model.find(self._sbml_tag('listOfReactions'))
        for reaction in self._reactions.iterfind(self._sbml_tag('reaction')):
            reaction_name = reaction.get('name')
            reaction_id = reaction.get('id') if self._level > 1 else reaction_name
            reaction_rev = reaction.get('reversible', 'true') == 'true'

            left = []
            for species_id, value in self._parse_species_references(reaction, 'listOfReactants'):
                if species_id not in self._model_compounds:
                    raise ParseError('Reaction {} references non-existent species {}'.format(reaction_id, species_id))
                species_name, species_comp = self._model_compounds[species_id]
                left.append((Compound(species_id, compartment=species_comp), value))

            right = []
            for species_id, value in self._parse_species_references(reaction, 'listOfProducts'):
                if species_id not in self._model_compounds:
                    raise ParseError('Reaction {} references non-existent species {}'.format(reaction_id, species_id))
                species_name, species_comp = self._model_compounds[species_id]
                right.append((Compound(species_id, compartment=species_comp), value))

            # Add reaction to database
            direction = Reaction.Bidir if reaction_rev else Reaction.Right
            self._database.set_reaction(reaction_id, Reaction(direction, left, right))

    def _parse_species_references(self, root, name):
        '''Yield species id and parsed value for a speciesReference list'''
        for species in root.iterfind('./{}/{}'.format(self._sbml_tag(name), self._sbml_tag('speciesReference'))):
            species_id = species.get('species')
            species_value = species.get('stoichiometry', 1)

            value = Decimal(species_value)
            if value % 1 == 0:
                value = int(value)

            yield species_id, value

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
        '''Name of compound'''
        if compound.name not in self._model_compounds:
            raise ValueError('Unknown compound: {}'.format(compound))
        name, comp = self._model_compounds[compound.name]
        return name

    def get_reaction_notes_elements(self):
        '''Yield tuples of reaction ids, and notes as ElementTree elements'''
        for reaction in self._reactions.iterfind(self._sbml_tag('reaction')):
            reaction_id = reaction.get('id')
            yield reaction_id, reaction.find(self._sbml_tag('notes'))

    def get_compound_notes_elements(self):
        '''Yield tuples of compound ids, and notes as ElementTree elements'''
        for compound in self._species.iterfind(self._sbml_tag('species')):
            compound_id = compound.get('id')
            yield compound_id, compound.find(self._sbml_tag('notes'))
