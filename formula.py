def parse(s):
    '''Parse a formula string (e.g. C6H10O2)

    >>> import pprint
    >>> pprint.pprint(parse_formula('H2O2'))
    {'H': 2, 'O': 2}
    >>> pprint.pprint(parse_formula('H2O'))
    {'H': 2, 'O': 1}
    >>> pprint.pprint(parse_formula('Zn'))
    {'Zn': 1}
    >>> pprint.pprint(parse_formula('C2H5NO2'))
    {'C': 2, 'H': 5, 'N': 1, 'O': 2}
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

    return formula

def multiply(formula, n):
    return { key: value*n for key, value in formula.items() }

def sum(f1, f2):
    '''Return sum of two formulas'''
    formula = {}
    for f in (f1, f2):
        for key, value in f.items():
            if key in formula:
                formula[key] += value
            else:
                formula[key] = value
    return formula

if __name__ == '__main__':

    import doctest
    doctest.testmod()
