
"""Module for reading and writing internal table-based formats

These formats are all space-separated text files. Empty lines are
ignored, just as comments starting with pound (#).
"""

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

def parse_limits_file(f):
    """Parse a space-separated file containing reaction flux limits

    The first column contains reaction IDs while the second column contains
    the lower flux limits. The third column is optional and contains the
    upper flux limit.
    """

    for line in f:
        line, _, comment = line.partition('#')
        line = line.strip()
        # TODO Comments can start with an asterisk to remain
        # compatible with GAMS files. Can be removed when
        # compatibility is no longer needed.
        if line == '' or line[0] == '*':
            continue

        # A line can specify lower limit only (useful for
        # exchange reactions), or both lower and upper limit.
        fields = line.split(None)
        if len(fields) == 2:
            reaction_id, lower = fields
            yield reaction_id, float(lower), None
        elif len(fields) == 3:
            reaction_id, lower, upper = fields
            yield reaction_id, float(lower), float(upper)
        else:
            raise ValueError('Malformed reaction limit: {}'.format(fields))

def parse_model_file(f):
    """Parse a file containing a list of reactions"""

    for line in f:
        line, _, comment = line.partition('#')
        line = line.strip()
        if line == '':
            continue
        yield line
