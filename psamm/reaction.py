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

"""Definitions related to reaction equations and parsing of such equations."""

import functools


@functools.total_ordering
class Compound(object):
    """Represents a compound in a reaction equation

    A compound is a named entity in the reaction equations representing a
    chemical compound. A compound can represent a generalized chemical entity
    (e.g. polyphosphate) and the arguments can be used to instantiate a
    specific chemical entity (e.g. polyphosphate(3)) by passing a number as an
    argument or a partially specified entity by passing an expression (e.g.
    polyphosphate(n)).
    """

    def __init__(self, name, compartment=None, arguments=()):
        self._name = str(name)
        self._compartment = None if compartment is None else str(compartment)
        self._arguments = tuple(arguments)

    @property
    def name(self):
        """Name of compound"""
        return self._name

    @property
    def compartment(self):
        """Compartment of compound"""
        return self._compartment

    @property
    def arguments(self):
        """Expression argument for generalized compounds"""
        return self._arguments

    def translate(self, func):
        """Translate compound name using given function

        >>> Compound('Pb').translate(lambda x: x.lower())
        Compound('pb')
        """
        return self.__class__(
            func(self._name), self._compartment, self._arguments)

    def in_compartment(self, compartment):
        """Return an instance of this compound in the specified compartment

        >>> Compound('H+').in_compartment('e')
        Compound('H+', 'e')
        """
        return self.__class__(self._name, compartment, self._arguments)

    def __eq__(self, other):
        return (isinstance(other, Compound) and
                self._name == other._name and
                self._compartment == other._compartment and
                self._arguments == other._arguments)

    def __ne__(self, other):
        return not self == other

    def __lt__(self, other):
        if isinstance(other, Compound):
            def tuple_repr(c):
                comp = c._compartment if c._compartment is not None else ''
                return (c._name, comp, c._arguments)
            return tuple_repr(self) < tuple_repr(other)
        return NotImplemented

    def __hash__(self):
        return (hash('Compound') ^ hash(self._name) ^
                hash(self._compartment) ^ hash(self._arguments))

    def __str__(self):
        """String representation of compound

        >>> str(Compound('Phosphate'))
        'Phosphate'
        >>> str(Compound('Phosphate', 'e'))
        'Phosphate[e]'
        >>> str(Compound('Polyphosphate', None, [Expression('n')]))
        'Polyphosphate(n)'
        >>> str(Compound('Polyphosphate', 'p', [Expression('n')]))
        'Polyphosphate(n)[p]'
        """
        s = self._name
        if len(self._arguments) > 0:
            s += '({})'.format(', '.join(str(a) for a in self._arguments))
        if self._compartment is not None:
            s += '[{}]'.format(self._compartment)
        return s

    def __repr__(self):
        def str_repr(*args):
            return 'Compound({})'.format(', '.join(repr(a) for a in args))

        if len(self._arguments) == 0:
            if self._compartment is None:
                return str_repr(self._name)
            return str_repr(self._name, self._compartment)
        return str_repr(self._name, self._compartment, self._arguments)


class Reaction(object):
    """Reaction equation representation

    Each compound is associated with a stoichiometric value.
    """

    Bidir = '<=>'
    Left = '<='
    Right = '=>'

    def __init__(self, direction, left, right):
        if direction not in (Reaction.Bidir, Reaction.Left, Reaction.Right):
            raise ValueError('Invalid direction: {}'.format(direction))
        self._direction = direction
        self._left = tuple(left)
        self._right = tuple(right)

    @property
    def direction(self):
        """Direction of reaction equation"""
        return self._direction

    @property
    def left(self):
        """Compounds on the left-hand side of the reaction equation."""
        return self._left

    @property
    def right(self):
        """Compounds on the right-hand side of the reaction equation."""
        return self._right

    @property
    def compounds(self):
        """Sequence of compounds on both sides of the reaction equation

        The sign of the stoichiometric values reflect whether the compound is
        on the left-hand side (negative) or the right-hand side (positive).
        """
        return tuple((c, -v) for c, v in self._left) + self._right

    def normalized(self):
        """Return normalized reaction

        The normalized reaction will have direction Bidir or Right.
        """

        if self._direction == Reaction.Left:
            direction = Reaction.Right
            left = self._right
            right = self._left
        else:
            direction = self._direction
            left = self._left
            right = self._right

        return Reaction(direction, left, right)

    def translated_compounds(self, translate):
        """Return reaction where compound names have been translated.

        For each compound the translate function is called with the compound
        name and the returned value is used as the new compound name. A new
        reaction is returned with the substituted compound names.
        """
        left = ((compound.translate(translate), count)
                for compound, count in self._left)
        right = ((compound.translate(translate), count)
                 for compound, count in self._right)

        return Reaction(self._direction, left, right)

    def copy(self):
        """Returns a distinct copy as a new Reaction object."""
        return self.__class__(self._direction, self._left, self._right)

    def __str__(self):
        # Use the same format as ModelSEED
        def format_compound(compound, count):
            """Format compound"""
            cpdspec = str(compound)
            if count != 1:
                return '({}) |{}|'.format(count, cpdspec)
            return '|{}|'.format(cpdspec)

        def format_compound_list(cmpds):
            """Format compound list"""
            return ' + '.join(format_compound(compound, count)
                              for compound, count in cmpds)

        return '{} {} {}'.format(
            format_compound_list(self._left),
            '?' if self._direction == '' else self._direction,
            format_compound_list(self._right))

    def __repr__(self):
        return 'Reaction({}, {}, {})'.format(
            repr(self._direction), repr(self._left), repr(self._right))

    def __eq__(self, other):
        """Indicate equality of self and other"""
        return (self._direction == other._direction and
                self._left == other._left and
                self._right == other._right)

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return (hash('Reaction') ^ hash(self._direction) ^
                hash(self._left) ^ hash(self._right))


if __name__ == '__main__':
    import doctest
    doctest.testmod()
