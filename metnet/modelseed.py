
'''Module related to loading ModelSEED database files'''

import csv
import re

class CompoundEntry(object):
    '''Representation of entry in a ModelSEED compound table'''

    def __init__(self, id, names, formula):
        self._id = id
        self._names = list(names)
        self._formula = formula

        # Find shortest name
        # Usually the best name is the shortest one,
        # but in the case of ties, the best one is
        # usually the last one to appear in the list.
        # Since the sort is stable we can obtain this
        # by first reversing the list, then sorting
        # by length.
        self._name = None
        if len(self._names) > 0:
            self._name = sorted(reversed(self._names), key=lambda x: len(x))[0]

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        return self._name

    @property
    def names(self):
        return iter(self._names)

    @property
    def formula(self):
        return self._formula


def parse_compound_file(f):
    '''Iterate over the compound entries in the given file'''

    f.readline() # Skip header
    for row in csv.reader(f, delimiter='\t'):
        compound_id, names, formula = row[:3]
        names = names.split(',<br>')

        # ModelSEED sometimes uses an asterisk and number at
        # the end of formulas. This seems to have a similar
        # meaning as '(...)n'.
        m = re.match(r'^(.*)\*(\d*)$', formula)
        if m is not None:
            if m.group(2) != '':
                formula = '({}){}'.format(m.group(1), m.group(2))
            else:
                formula = '({})n'.format(m.group(1))

        formula = formula.strip()
        if formula == '' or formula == 'noformula':
            formula = None

        yield CompoundEntry(compound_id, names, formula)