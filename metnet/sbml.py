
'''Parser for SBML model files'''

import xml.etree.ElementTree as ET
from decimal import Decimal

from .metabolicmodel import DictDatabase
from .reaction import Reaction, Compound

SBML_NS = 'http://www.sbml.org/sbml/level2'

def sbml_name(tag):
    '''Prepend namespace to SBML tag name'''
    return '{{{}}}{}'.format(SBML_NS, tag)

def parse_species_references(root, tag_name):
    '''Yield species id and parsed value for a speciesReference list'''
    for species in root.findall('./{}/{}'.format(tag_name, sbml_name('speciesReference'))):
        species_id = species.get('species')
        species_value = species.get('stoichiometry', 1)

        value = Decimal(species_value)
        if value % 1 == 0:
            value = int(value)

        yield species_id, value

def parse_sbml_file(file):
    '''Parse SBML file and return database of reactions'''

    # Parse XML file
    tree = ET.parse(file)
    root = tree.getroot()

    # Only one model is allowed in each file
    model = root.find('./{}'.format(sbml_name('model')))

    # Compartments
    model_compartments = {}
    for compartment in model.findall('./{}/{}'.format(sbml_name('listOfCompartments'), sbml_name('compartment'))):
        comp_id = compartment.get('id')
        comp_name = compartment.get('name')
        model_compartments[comp_id] = comp_name

    # Compounds
    model_compounds = {}
    for species in model.findall('./{}/{}'.format(sbml_name('listOfSpecies'), sbml_name('species'))):
        species_id = species.get('id')
        species_name = species.get('name')
        species_comp = species.get('compartment')
        model_compounds[species_id] = species_name, species_comp

    database = DictDatabase()

    # Reactions
    for reaction in model.findall('./{}/{}'.format(sbml_name('listOfReactions'), sbml_name('reaction'))):
        reaction_id = reaction.get('id')
        reaction_name = reaction.get('name')
        reaction_rev = reaction.get('reversible', 'true').lower() in ('true', 'yes', '1')

        left = []
        for species_id, value in parse_species_references(reaction, sbml_name('listOfReactants')):
            species_name, species_comp = model_compounds[species_id]
            left.append((Compound(species_id, compartment=species_comp), value))

        right = []
        for species_id, value in parse_species_references(reaction, sbml_name('listOfProducts')):
            species_name, species_comp = model_compounds[species_id]
            right.append((Compound(species_id, compartment=species_comp), value))

        # Add reaction to database
        direction = Reaction.Bidir if reaction_rev else Reaction.Right
        database.set_reaction(reaction_id, Reaction(direction, left, right))

    return database