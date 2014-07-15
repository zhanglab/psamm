
import re

class ParseError(Exception):
    '''Exception used to signal errors while parsing'''
    pass

class CompoundEntry(object):
    '''Representation of entry in KEGG compound file'''

    def __init__(self, values):
        self.values = dict(values)
        if 'entry' not in values:
            raise ParseError('Missing compound identifier')
        self._id, _ = values['entry'][0].split(None, 1)

    @property
    def id(self):
        return self._id

    @property
    def names(self):
        if 'name' in self.values:
            for line in self.values['name']:
                for name in line.rstrip(';').split(';'):
                    yield name.strip()

    @property
    def reactions(self):
        if 'reaction' in self.values:
            for line in self.values['reaction']:
                for rxnid in line.split():
                    yield rxnid

    @property
    def enzymes(self):
        if 'enzyme' in self.values:
            for line in self.values['enzyme']:
                for enzyme in line.split():
                    yield enzyme

    @property
    def formula(self):
        if 'formula' not in self.values:
            return None
        return self.values['formula'][0]

    @property
    def exact_mass(self):
        if 'exact_mass' not in self.values:
            return None
        return float(self.values['exact_mass'][0])

    @property
    def mol_weight(self):
        if 'mol_weight' not in self.values:
            return None
        return float(self.values['mol_weight'][0])

    @property
    def pathways(self):
        if 'pathway' in self.values:
            for line in self.values['pathway']:
                pathway, name = line.split(None, 1)
                yield pathway, name

    @property
    def dblinks(self):
        if 'dblinks' in self.values:
            for line in self.values['dblinks']:
                database, entry = line.split(':', 1)
                yield database.strip(), entry.strip()

    def __getitem__(self, name):
        if name not in self.values:
            raise AttributeError('Attribute does not exist: {}'.format(name))
        return self._values[name]

    def __repr__(self):
        return '<CompoundEntry "{}">'.format(self.id)

def parse_compound_file(f):
    '''Iterate over the compound entries in the given file'''

    section_id = None
    compound = {}
    for line in f:
        if line.strip() == '///':
            # End of compound
            yield CompoundEntry(compound)
            compound = {}
            section_id = None
        else:
            # Look for beginning of section
            m = re.match(r'([A-Z_]+)\s+(.*)', line.rstrip())
            if m is not None:
                section_id = m.group(1).lower()
                compound[section_id] = [m.group(2)]
            elif section_id is not None:
                compound[section_id].append(line.strip())
            else:
                raise ParseError('Missing section identifier')
