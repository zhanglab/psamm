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

from __future__ import unicode_literals

import re

import six
from six import text_type, string_types

from ..util import FrozenOrderedSet


@six.python_2_unicode_compatible
class Variable(object):
    """Represents a variable in a boolean expression"""

    def __init__(self, symbol):
        self._symbol = text_type(symbol)

    @property
    def symbol(self):
        return self._symbol

    def __repr__(self):
        return str('Variable({})').format(repr(self._symbol))

    def __str__(self):
        return self._symbol

    def __eq__(self, other):
        return isinstance(other, Variable) and self._symbol == other._symbol

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash('Variable') ^ hash(self._symbol)


class _OperatorTerm(object):
    """Composite operator term."""
    def __init__(self, *args):
        terms = list()
        for arg in args:
            if isinstance(arg, self.__class__):
                terms.extend(arg)
            elif isinstance(arg, (bool, Variable, _OperatorTerm)):
                terms.append(arg)
            else:
                raise ValueError('Invalid term: {!r}'.format(arg))
        self._terms = FrozenOrderedSet(terms)

    def __iter__(self):
        return iter(self._terms)

    def __hash__(self):
        return hash(self._terms)

    def __eq__(self, other):
        return (isinstance(other, self.__class__) and
                self._terms == other._terms)

    def __ne__(self, other):
        return not self == other


class And(_OperatorTerm):
    """Represents an AND term in an expression."""


class Or(_OperatorTerm):
    """Represents an OR term in an expression."""


class SubstitutionError(Exception):
    """Error substituting into expression."""


@six.python_2_unicode_compatible
class Expression(object):
    """Boolean expression representation.

    The expression can be constructed from an expression string of
    variables, operators ("and", "or") and parenthesis groups. For example,

    >>> e = Expression('a and (b or c)')
    """

    def __init__(self, arg, _vars=None):
        if isinstance(arg, (_OperatorTerm, Variable, bool)):
            self._root = arg
        elif isinstance(arg, string_types):
            self._root = _parse_expression(arg)
        else:
            raise TypeError('Unexpected arguments to Expression: {}'.format(
                repr(arg)))

        # If present use _vars to create the set of variables directly. This
        # saves a loop over the tree nodes when the variables are already
        # known.
        if _vars is None:
            variables = []
            if isinstance(self._root, (_OperatorTerm, Variable)):
                stack = [self._root]
                while len(stack) > 0:
                    term = stack.pop()
                    if isinstance(term, Variable):
                        variables.append(term)
                    elif isinstance(term, _OperatorTerm):
                        stack.extend(reversed(list(term)))
                    else:
                        raise ValueError(
                            'Invalid node in expression tree: {!r}'.format(
                                term))

            self._variables = FrozenOrderedSet(variables)
        else:
            self._variables = FrozenOrderedSet(_vars)

    @property
    def root(self):
        """Return root term, variable or boolean of the expression."""
        return self._root

    @property
    def variables(self):
        """Immutable set of variables in the expression."""
        return self._variables

    def has_value(self):
        """Return True if the expression has no variables."""
        return isinstance(self._root, bool)

    @property
    def value(self):
        """The value of the expression if fully evaluated."""
        if not self.has_value():
            raise ValueError('Expression is not fully evaluated')
        return self._root

    def substitute(self, mapping):
        """Substitute variables using mapping function."""
        next_terms = iter([self._root])
        output_stack = []
        current_type = None
        terms = []
        variables = []
        term = None

        while True:
            try:
                term = next(next_terms)
            except StopIteration:
                term = None

            if term is None:
                if current_type is None:
                    term = terms[0]
                    break
                else:
                    if len(terms) == 0:
                        if current_type == And:
                            term = True
                        elif current_type == Or:
                            term = False
                    elif len(terms) == 1:
                        term = terms[0]
                    else:
                        term = current_type(*terms)
                current_type, next_terms, terms = output_stack.pop()
            else:
                if isinstance(term, _OperatorTerm):
                    output_stack.append((current_type, next_terms, terms))
                    current_type = term.__class__
                    terms = []
                    next_terms = iter(term)
                    continue

            # Substitute variable
            if isinstance(term, Variable):
                term = mapping(term)
                if not isinstance(term, (_OperatorTerm, Variable, bool)):
                    raise SubstitutionError(
                        'Expected Variable or bool from substitution,'
                        ' got: {!r}'.format(term))

            # Check again after substitution
            if isinstance(term, Variable):
                variables.append(term)

            # Short circuit with booleans
            while isinstance(term, bool):
                if current_type == And:
                    if not term:
                        current_type, next_terms, terms = output_stack.pop()
                        continue
                    else:
                        break
                elif current_type == Or:
                    if term:
                        current_type, next_terms, terms = output_stack.pop()
                        continue
                    else:
                        break
                else:
                    terms.append(term)
                    break
            else:
                terms.append(term)

        return self.__class__(term, _vars=variables)

    def __repr__(self):
        if isinstance(self._root, bool):
            arg = self._root
        else:
            arg = text_type(self)
        return str('{}({!r})').format(self.__class__.__name__, arg)

    def __str__(self):
        next_terms = iter([self._root])
        output_stack = []
        current_type = None
        terms = []
        term = None

        while True:
            try:
                term = next(next_terms)
            except StopIteration:
                term = None

            if term is None:
                if current_type is None:
                    term = terms[0]
                    break
                elif current_type == And:
                    term = ' and '.join(t for t in terms)
                elif current_type == Or:
                    term = ' or '.join(t for t in terms)
                current_type, next_terms, terms = output_stack.pop()

                # Break on None here to avoid wrapping the outermost term in
                # parentheses.
                if current_type is None:
                    break

                terms.append('(' + term + ')')
            else:
                if isinstance(term, _OperatorTerm):
                    output_stack.append((current_type, next_terms, terms))
                    current_type = term.__class__
                    terms = []
                    next_terms = iter(term)
                else:
                    terms.append(text_type(term))

        return term

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self._root == other._root
        return NotImplemented

    def __ne__(self, other):
        if isinstance(other, self.__class__):
            return not self == other
        return NotImplemented


class ParseError(Exception):
    """Signals error parsing boolean expression."""

    def __init__(self, *args, **kwargs):
        self._span = kwargs.pop('span', None)
        super(ParseError, self).__init__(*args, **kwargs)

    @property
    def indicator(self):
        if self._span is None:
            return None
        pre = ' ' * self._span[0]
        ind = '^' * max(1, self._span[1] - self._span[0])
        return pre + ind


def _parse_expression(s):
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
            raise ParseError('Invalid token in expression string: {}'.format(
                repr(match.group(0))), span=(match.start(), match.end()))
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
                raise ParseError(
                    'Unbalanced parenthesis group in expression',
                    span=(match.start(), match.end()))
            if group_pairs[group_end] != clause_symbol:
                raise ParseError(
                    'Group started with {} ended with {}'.format(
                        clause_symbol, group_end),
                    span=(match.start(), match.end()))
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
            raise ParseError(
                'Invalid token in expression string: {!r}'.format(
                    match.group(0)),
                span=(match.start(), match.end()))

    if len(clause_stack) > 0:
        raise ParseError('Unbalanced parenthesis group in expression')

    expr = operators[clause_operator](*current_clause)
    return expr
