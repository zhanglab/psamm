
"""Module for reading and writing internal table-based formats"""

from ..reaction import Reaction, Compound, ModelSEED

def parse_reaction_file(f):
    """Parse a space-separated file containing reaction IDs and definitions

    The reaction definitions are parsed as ModelSEED format.
    """

    for line in f:
        line, _, comment = line.partition('#')
        line = line.strip()
        if line == '':
            continue
        reaction_id, equation = line.split(None, 1)
        reaction = ModelSEED.parse(equation).normalized()
        yield reaction_id, reaction
