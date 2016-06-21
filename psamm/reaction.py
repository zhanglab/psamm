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

from __future__ import unicode_literals

import functools
import enum
import numbers
from collections import Counter

import six
from six import text_type, iteritems


@six.python_2_unicode_compatible
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
        self._name = text_type(name)
        self._compartment = None
        if compartment is not None:
            self._compartment = text_type(compartment)
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
            s += '({})'.format(
                ', '.join(text_type(a) for a in self._arguments))
        if self._compartment is not None:
            s += '[{}]'.format(self._compartment)
        return s

    def __repr__(self):
        def str_repr(*args):
            return str('Compound({})').format(', '.join(repr(a) for a in args))

        if len(self._arguments) == 0:
            if self._compartment is None:
                return str_repr(self._name)
            return str_repr(self._name, self._compartment)
        return str_repr(self._name, self._compartment, self._arguments)


class Direction(enum.Enum):
    """Directionality of reaction equation."""
    Forward = False, True
    Reverse = True, False
    Both = True, True

    Right = False, True
    Left = True, False
    Bidir = True, True

    @property
    def forward(self):
        """Whether this direction includes forward direction."""
        return self.value[1]

    @property
    def reverse(self):
        """Whether this direction includes reverse direction."""
        return self.value[0]

    def flipped(self):
        """Return the flipped version of this direction."""
        forward, reverse = self.value
        return self.__class__((reverse, forward))

    @property
    def symbol(self):
        """Return string symbol for direction."""
        if self == Direction.Forward:
            return '=>'
        elif self == Direction.Reverse:
            return '<='
        else:
            return '<=>'


@six.python_2_unicode_compatible
class Reaction(object):
    """Reaction equation representation.

    Each compound is associated with a stoichiometric value and the reaction
    has a :class:`.Direction`. The reaction is created in one of the three
    following ways.

    It can be created from a direction and two iterables of compound,
    value pairs representing the left-hand side and the right-hand side of
    the reaction:

    >>> r = Reaction(Direction.Both, [(Compound('A'), 1), (Compound('B', 2))],
                                     [(Compound('C'), 1)])
    >>> str(r)
    '|A| + (2) |B| <=> |C|'

    It can also be created from a single dict or iterable of compound, value
    pairs where the left-hand side compounds have negative values and the
    right-hand side compounds have positive values:

    >>> r = Reaction(Direction.Forward, {
            Compound('A'): -1,
            Compound('B'): -2,
            Compound('C'): 1
    })
    >>> str(r)
    '|A| + (2) |B| <=> |C|'

    Lastly, the reaction can be created from an existing reaction object,
    creating a copy of that reaction.

    >>> r = Reaction(Direction.Forward, {Compound('A'): -1, Compound('B'): 1})
    >>> r2 = Reaction(r)
    >>> str(r2)
    '|A| => |B|'

    Reactions can be added to produce combined reactions.

    >>> r = Reaction(Direction.Forward, {Compound('A'): -1, Compound('B'): 1})
    >>> s = Reaction(Direction.Forward, {Compound('B'): -1, Compound('C'): 1})
    >>> str(r + s)
    '|A| => |C|'

    Reactions can also be multiplied by a number to produce a new reaction
    with scaled stoichiometry.

    >>> r = Reaction(Direction.Forward, {Compound('A'): -1, Compound('B'): 2})
    >>> str(2 * r)
    '(2) |A| => (4) |B|'

    Multiplying with a negative value will also flip the reaction, and as a
    special case, negating a reaction will simply flip it.

    >>> r = Reaction(Direction.Forward, {Compound('A'): -1, Compound('B'): 2})
    >>> str(r)
    '|A| => (2) |B|'
    >>> str(-r)
    '(2) |B| <= |A|'
    """

    def __init__(self, *args):
        if len(args) == 1:
            # Initialize from Reaction object
            if not isinstance(args[0], self.__class__):
                raise TypeError('Single argument must be of type {}'.format(
                    self.__class__.__name__))

            self._direction = args[0].direction
            self._left = tuple(args[0].left)
            self._right = tuple(args[0].right)
        elif len(args) == 2:
            # Initialize from direction and dict or single iterable
            if not isinstance(args[0], Direction):
                raise TypeError('First argument must be a Direction')
            self._direction = args[0]

            values = args[1]
            if isinstance(values, dict):
                values = iteritems(values)

            left, right = [], []
            for compound, value in values:
                if not isinstance(value, numbers.Number):
                    raise TypeError('Values must be numeric')
                if value < 0:
                    left.append((compound, -value))
                elif value > 0:
                    right.append((compound, value))

            self._left = tuple(left)
            self._right = tuple(right)
        elif len(args) == 3:
            # Initialize from direction and two iterables
            if not isinstance(args[0], Direction):
                raise TypeError('First argument must be a Direction')
            self._direction = args[0]

            left, right = [], []
            for compound, value in args[1]:
                if isinstance(value, numbers.Number) and value < 0:
                    raise ValueError('Value must not be negative')
                elif value != 0:
                    left.append((compound, value))

            for compound, value in args[2]:
                if isinstance(value, numbers.Number) and value < 0:
                    raise ValueError('Value must not be negative')
                elif value != 0:
                    right.append((compound, value))

            self._left = tuple(left)
            self._right = tuple(right)
        else:
            raise TypeError('Too many arguments (one, two or three'
                            ' arguments required)')

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

        The normalized reaction will be bidirectional or a forward reaction
        (i.e. reverse reactions are flipped).
        """

        if self._direction == Direction.Reverse:
            return self.__class__(
                self._direction.flipped(), self._right, self._left)

        return self

    def translated_compounds(self, translate):
        """Return reaction where compound names have been translated.

        For each compound the translate function is called with the compound
        name and the returned value is used as the new compound name. A new
        reaction is returned with the substituted compound names.
        """
        compounds = ((compound.translate(translate), value)
                     for compound, value in self.compounds)
        return Reaction(self._direction, compounds)

    def __str__(self):
        # Use the same format as ModelSEED
        def format_compound(compound, count):
            """Format compound"""
            cpdspec = text_type(compound)
            if count != 1:
                return '({}) |{}|'.format(count, cpdspec)
            return '|{}|'.format(cpdspec)

        def format_compound_list(cmpds):
            """Format compound list"""
            return ' + '.join(format_compound(compound, count)
                              for compound, count in cmpds)

        return '{} {} {}'.format(
            format_compound_list(self._left),
            self._direction.symbol,
            format_compound_list(self._right))

    def __repr__(self):
        return str('Reaction({}, {}, {})').format(
            repr(self._direction), repr(self._left), repr(self._right))

    def __eq__(self, other):
        """Indicate equality of self and other"""
        return (isinstance(other, self.__class__) and
                self._direction == other._direction and
                self._left == other._left and
                self._right == other._right)

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return (hash('Reaction') ^ hash(self._direction) ^
                hash(self._left) ^ hash(self._right))

    def __add__(self, other):
        if isinstance(other, Reaction):
            reverse = self.direction.reverse and other.direction.reverse
            forward = self.direction.forward and other.direction.forward
            if not reverse and not forward:
                raise ValueError('Reactions have incompatible directions')

            direction = Direction((reverse, forward))
            values = Counter(dict(self.compounds))
            values.update(dict(other.compounds))
            return Reaction(direction, values)

        return NotImplemented

    def __sub__(self, other):
        return self + -other

    def __neg__(self):
        return -1 * self

    def __mul__(self, other):
        if isinstance(other, numbers.Number):
            direction = (
                self.direction if other > 0 else self.direction.flipped())
            return self.__class__(
                direction, ((c, other * v) for c, v in self.compounds))

        return NotImplemented

    def __rmul__(self, other):
        return self * other
