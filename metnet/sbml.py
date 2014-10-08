
'''Parser for SBML model files'''

import xml.etree.ElementTree as ET
from decimal import Decimal

from .metabolicmodel import MetabolicDatabase, DictDatabase
from .reaction import Reaction, Compound

SBML_NS = 'http://www.sbml.org/sbml/level2'

def sbml_name(tag):
    '''Prepend namespace to SBML tag name'''
    return '{{{}}}{}'.format(SBML_NS, tag)

def parse_species_references(root, tag_name):
    '''Yield species id and parsed value for a speciesReference list'''
    for species in root.iterfind('./{}/{}'.format(tag_name, sbml_name('speciesReference'))):
        species_id = species.get('species')
        species_value = species.get('stoichiometry', 1)

        value = Decimal(species_value)
        if value % 1 == 0:
            value = int(value)

        yield species_id, value

class SBMLDatabase(MetabolicDatabase):
    def __init__(self, file):
        # Parse SBML file
        tree = ET.parse(file)
        root = tree.getroot()
        self._model = root.find(sbml_name('model'))

        self._database = DictDatabase()

        # Compounds
        self._model_compounds = {}
        self._species = self._model.find(sbml_name('listOfSpecies'))
        for species in self._species.iterfind(sbml_name('species')):
            species_id = species.get('id')
            species_name = species.get('name')
            species_comp = species.get('compartment')
            self._model_compounds[species_id] = species_name, species_comp

        # Reactions
        self._reactions = self._model.find(sbml_name('listOfReactions'))
        for reaction in self._reactions.iterfind(sbml_name('reaction')):
            reaction_id = reaction.get('id')
            reaction_name = reaction.get('name')
            reaction_rev = reaction.get('reversible', 'true') == 'true'

            left = []
            for species_id, value in parse_species_references(reaction, sbml_name('listOfReactants')):
                species_name, species_comp = self._model_compounds[species_id]
                left.append((Compound(species_id, compartment=species_comp), value))

            right = []
            for species_id, value in parse_species_references(reaction, sbml_name('listOfProducts')):
                species_name, species_comp = self._model_compounds[species_id]
                right.append((Compound(species_id, compartment=species_comp), value))

            # Add reaction to database
            direction = Reaction.Bidir if reaction_rev else Reaction.Right
            self._database.set_reaction(reaction_id, Reaction(direction, left, right))

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
        for reaction in self._reactions.iterfind(sbml_name('reaction')):
            reaction_id = reaction.get('id')
            yield reaction_id, reaction.find(sbml_name('notes'))

    def get_compound_notes_elements(self):
        '''Yield tuples of compound ids, and notes as ElementTree elements'''
        for compound in self._species.iterfind(sbml_name('species')):
            compound_id = compound.get('id')
            yield compound_id, compound.find(sbml_name('notes'))
