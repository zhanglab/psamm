# This file is part of PSAMM.
#
# PSAMM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PSAMM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PSAMM.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2014-2015  Jon Lund Steffensen <jon_steffensen@uri.edu>

"""Representations of boolean expressions and variables.

These classes can be used to represent simple boolean
expressions and do evaluation with substitutions of
particular variables.
"""

import re

from six import text_type


class Variable(object):
    """Represents a variable in a boolean expression"""

    def __init__(self, symbol):
        self._symbol = text_type(symbol)

    @property
    def symbol(self):
        return self._symbol

    def substitute(self, mapping):
        """Substitute variables using mapping function"""
        return mapping(self)

    def variables(self):
        """Iterate over variables."""
        yield self

    def __repr__(self):
        return 'Variable({})'.format(repr(self._symbol))

    def __str__(self):
        return self._symbol

    def __eq__(self, other):
        return isinstance(other, Variable) and self._symbol == other._symbol

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash('Variable') ^ hash(self._symbol)


class And(object):
    """Represents a conjuction of boolean terms"""

    def __init__(self, *args):
        terms = set()
        for t in args:
            if isinstance(t, And):
                terms.update(t.terms)
            else:
                terms.add(t)
        self._terms = frozenset(terms)

    @property
    def terms(self):
        """Iterate over terms"""
        return iter(self._terms)

    def variables(self):
        """Iterate over variables."""
        for t in self._terms:
            if isinstance(t, Variable):
                yield t
            else:
                for u in t.variables():
                    yield u

    def substitute(self, mapping):
        """Substitute variables using mapping function"""
        result = []
        for t in self._terms:
            value = t.substitute(mapping)
            if isinstance(value, bool):
                if not value:
                    return False
            else:
                result.append(value)
        if len(result) == 0:
            return True
        elif len(result) == 1:
            return result[0]
        return And(*result)

    def __repr__(self):
        return 'Expression({})'.format(repr(text_type(self)))

    def __str__(self):
        def format_term(t):
            if isinstance(t, Variable):
                return text_type(t)
            return '({})'.format(t)
        return ' and '.join(format_term(t) for t in self._terms)

    def __eq__(self, other):
        if isinstance(other, Variable) and len(self._terms) == 1:
            return self._terms[0] == other
        return isinstance(other, And) and self._terms == other._terms

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash('And') ^ hash(self._terms)


class Or(object):
    """Represents a disjuction of boolean terms"""

    def __init__(self, *args):
        terms = set()
        for t in args:
            if isinstance(t, Or):
                terms.update(t.terms)
            else:
                terms.add(t)
        self._terms = frozenset(terms)

    @property
    def terms(self):
        """Iterate over terms"""
        return iter(self._terms)

    def variables(self):
        """Iterate over variables."""
        for t in self._terms:
            if isinstance(t, Variable):
                yield t
            else:
                for u in t.variables():
                    yield u

    def substitute(self, mapping):
        """Substitute variables using mapping function"""
        result = []
        for t in self._terms:
            value = t.substitute(mapping)
            if isinstance(value, bool):
                if value:
                    return True
            else:
                result.append(value)
        if len(result) == 0:
            return False
        elif len(result) == 1:
            return result[0]
        return Or(*result)

    def __repr__(self):
        return 'Expression({})'.format(repr(text_type(self)))

    def __str__(self):
        def format_term(t):
            if isinstance(t, Variable):
                return text_type(t)
            return '({})'.format(t)
        return ' or '.join(format_term(t) for t in self._terms)

    def __eq__(self, other):
        if isinstance(other, Variable) and len(self._terms) == 1:
            return self._terms[0] == other
        return isinstance(other, Or) and self._terms == other._terms

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash('Or') ^ hash(self._terms)


def Expression(s):  # noqa
    """Parse boolean expression containing and/or operators"""

    # Converters for opeartor clauses
    operators = {
        'and': And,
        'or': Or,
        None: lambda *args: args[0]
    }

    # Pairing of end group symbols with start group symbols
    group_pairs = {
        ')': '(',
        ']': '['
    }

    scanner = re.compile(r'''
        (\s+) |              # space
        (\(|\[) |            # group_start
        (\)|\]) |            # group_end
        ((?:or|and)\b) |     # operator
        ([^\s\(\)\[\]]+) |   # variable
        (\Z) |               # end
        (.)                  # error
    ''', re.DOTALL | re.VERBOSE | re.UNICODE | re.IGNORECASE)

    # Parsed using two states and a stack of open clauses
    # At state 0 (not expect_operator): Expect variable, or parenthesis group
    #  start.
    # At state 1 (expect_operator): Expect operator, parenthesis group end, or
    #  end.
    expect_operator = False
    clause_stack = []
    current_clause = []
    clause_operator = None
    clause_symbol = None

    def close():
        prev_op, prev_symbol, prev_clause = clause_stack.pop()
        prev_clause.append(operators[clause_operator](*current_clause))
        return prev_op, prev_symbol, prev_clause

    for match in re.finditer(scanner, s):
        (space, group_start, group_end, operator, variable, end,
            error) = match.groups()

        if error is not None:
            raise ValueError('Invalid token in expression string: {}'.format(
                repr(match.group(0))))
        elif space is not None:
            continue
        elif expect_operator and operator is not None:
            operator = operator.lower()
            if operator == 'and' and clause_operator != 'and':
                prev_term = current_clause.pop()
                clause_stack.append(
                    (clause_operator, clause_symbol, current_clause))
                current_clause = [prev_term]
            elif operator == 'or' and clause_operator == 'and':
                clause_operator, clause_symbol, current_clause = close()
            clause_operator = operator
            expect_operator = False
        elif expect_operator and group_end is not None:
            if clause_operator == 'and':
                clause_operator, clause_symbol, current_clause = close()
            if len(clause_stack) == 0:
                raise ValueError('Unbalanced parenthesis group in expression')
            if group_pairs[group_end] != clause_symbol:
                raise ValueError(
                    'Group started with {} ended with {}'.format(
                        clause_symbol, group_end))
            clause_operator, clause_symbol, current_clause = close()
        elif expect_operator and end is not None:
            if clause_operator == 'and':
                clause_operator, clause_symbol, current_clause = close()
        elif not expect_operator and variable is not None:
            current_clause.append(Variable(variable))
            expect_operator = True
        elif not expect_operator and group_start is not None:
            clause_stack.append(
                (clause_operator, clause_symbol, current_clause))
            current_clause = []
            clause_operator = None
            clause_symbol = group_start
        else:
            raise ValueError('Invalid token in expression string: {}'.format(
                repr(match.group(0))))

    if len(clause_stack) > 0:
        raise ValueError('Unbalanced parenthesis group in expression')

    expr = operators[clause_operator](*current_clause)
    return expr
