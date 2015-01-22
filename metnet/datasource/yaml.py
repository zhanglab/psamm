
"""Module for reading YAML based formats"""

from __future__ import absolute_import

import yaml

from ..reaction import Reaction, Compound
from . import modelseed

class ParseError(Exception):
    """Exception used to signal errors while parsing"""

def parse_reaction_list(reactions):
    """Parse a list of reactions obtained from a YAML file

    Yields tuples of reaction ID and reaction object."""

    for reaction_def in reactions:
        reaction_id = reaction_def.get('id')
        if reaction_id is None:
            raise ParseError('Reaction ID missing')

        def parse_compound_list(l):
            """Parse a list of reactants or metabolites"""
            for compound_def in l:
                compound_id = compound_def.get('id')
                if compound_id is None:
                    raise ParseError('Compound ID missing')

                value = compound_def.get('value')
                if value is None:
                    raise ParseError('Missing value for compound {}'.format(compound_id))

                compound_compartment = compound_def.get('compartment')
                if compound_compartment is None:
                    compound_compartment = compartment

                yield Compound(compound_id, compartment=compound_compartment), value

        if 'equation' in reaction_def:
            if any(key in reaction_def for key in ('compartment', 'reversible',
                                                   'left', 'right')):
                raise ParseError('Reaction contains ambiguous fields')
            equation = reaction_def.get('equation')
            reaction = modelseed.parse_reaction(equation).normalized()
        else:
            compartment = reaction_def.get('compartment', None)
            reversible = bool(reaction_def.get('reversible', True))
            left = reaction_def.get('left', [])
            right = reaction_def.get('right', [])
            if len(left) == 0 and len(right) == 0:
                raise ParseError('Reaction values are missing')

            reaction = Reaction(Reaction.Bidir if reversible else Reaction.Right,
                                parse_compound_list(left),
                                parse_compound_list(right))

        yield reaction_id, reaction

def parse_reaction_file(f):
    """Parse a file as a YAML-format list of reactions"""
    return parse_reaction_list(yaml.load(f))
