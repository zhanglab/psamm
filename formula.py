
'''Parser for chemical formulas'''

import pprint
from collections import defaultdict

class Formula(object):
    def __init__(self, values={}):
        self._values = dict(values)

    def __repr__(self):
        return 'Formula({})'.format(pprint.pformat(self._values))

    def __add__(self, other):
        '''Sum of self and other

        >>> Formula({'H': 2, 'O': 1}) + Formula({'N': 1, 'O': 2})
        Formula({'H': 2, 'N': 1, 'O': 3})
        '''
        result = defaultdict(int)
        for f in (self._values, other._values):
            for key, value in f.items():
                result[key] += value
        return Formula(result)

    def __mul__(self, other):
        '''Multiply each count by other

        >>> Formula({'H': 2, 'O': 1}) * 4
        Formula({'H': 8, 'O': 4})'''

        return Formula({ key: value*other for key, value in self._values.items() })

    def __rmul__(self, other):
        '''Multiply each count by other

        >>> 2 * Formula({'H': 2, 'O': 1})
        Formula({'H': 4, 'O': 2})'''

        return self * other

    def __eq__(self, other):
        '''Equality of self and other

        >>> Formula({'H': 2, 'O': 1}) == Formula({'O': 1, 'H': 2})
        True
        >>> Formula({'Au': 1}) == Formula({'Ag': 1})
        False'''

        return self._values == other._values

    def __ne__(self, other):
        '''Not equality of self and other

        >>> Formula({'Au': 1}) != Formula({'Ag': 1})
        True'''

        return self._values != other._values

    @classmethod
    def parse(cls, s):
        '''Parse a formula string (e.g. C6H10O2)

        >>> Formula.parse('H2O2')
        Formula({'H': 2, 'O': 2})
        >>> Formula.parse('H2O')
        Formula({'H': 2, 'O': 1})
        >>> Formula.parse('Zn')
        Formula({'Zn': 1})
        >>> Formula.parse('C2H5NO2')
        Formula({'C': 2, 'H': 5, 'N': 1, 'O': 2})
        '''

        formula = {}
        looking_for_number = True
        has_read_uppercase = False

        number = ''
        atom = ''

        for i in range(len(s)):
            if looking_for_number:
                if s[i].isdigit():
                    looking_for_number = False
                    number += s[i]
                elif s[i].islower():
                    atom += s[i]
                else:
                    if has_read_uppercase:
                        # Implicit 1
                        formula[atom] = 1
                        atom = s[i]
                        number = ''
                    else:
                        atom += s[i]
                        has_read_uppercase = True
            else:
                if s[i].isdigit():
                    number += s[i]
                else:
                    formula[atom] = int(number)
                    atom = s[i]
                    number = ''
                    looking_for_number = True
                    has_read_uppercase = True

        formula[atom] = int(number) if number != '' else 1
        return Formula(formula)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
