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

"""Reaction parser."""

import re
from decimal import Decimal
import enum

from six import raise_from

from psamm.reaction import Reaction, Compound, Direction
from psamm.expression import affine


class ParseError(Exception):
    """Error raised when parsing reaction fails."""

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


@enum.unique
class _ReactionToken(enum.Enum):
    Space = 0
    Quoted = 1
    Group = 2
    Plus = 3
    Arrow = 4
    Other = 5
    End = 6


class ReactionParser(object):
    """Parser of reactions equations.

    This parser recognizes:

    * Global compartment specification as a prefix (when ``parse_global`` is
      ``True``)

    * Configurable reaction arrow tokens (``arrows``)

    * Compounds quoted by pipe (``|``) (required only if the compound name
      includes a space)

    * Compound counts that are affine expressions.
    """

    def __init__(self, arrows=None, parse_global=False):
        if arrows is None:
            arrows = (
                ('<=>', Direction.Both),
                ('=>', Direction.Forward),
                ('<=', Direction.Reverse)
            )

        arrow_p = '|'.join(re.escape(arrow) for arrow, _ in arrows)
        self._arrows = dict(arrows)
        self._tokens = (
            (r'\s+', _ReactionToken.Space),
            (r'\|[^|]+\|', _ReactionToken.Quoted),
            (r'\([^)]+\)(?=\Z|\s|\|)', _ReactionToken.Group),
            (r'\+(?=\s)', _ReactionToken.Plus),
            (arrow_p + r'(?=\Z|\s)', _ReactionToken.Arrow),
            (r'\S+', _ReactionToken.Other),
            (r'\Z', _ReactionToken.End),
            (r'.', None)
        )

        self._scanner = re.compile(
            '|'.join('(' + pattern + ')' for pattern, _ in self._tokens),
            re.DOTALL | re.UNICODE)

        self._parse_global = parse_global

    def parse(self, s):
        """Parse reaction string."""

        global_comp = None
        if self._parse_global:
            # Split by colon for global compartment information
            m = re.match(r'^\s*\[(\w+)\]\s*:\s*(.*)', s)
            if m is not None:
                global_comp = m.group(1)
                s = m.group(2)

        expect_operator = False
        direction = None
        left = []
        right = []
        current_side = left
        saved_token = None

        def tokenize():
            for match in re.finditer(self._scanner, s):
                for i, group in enumerate(match.groups()):
                    if group is not None:
                        token = self._tokens[i][1]
                        span = match.start(), match.end()
                        if token is None:
                            raise ParseError(
                                'Invalid token in expression string:'
                                ' {!r}'.format(group), span=span)
                        yield self._tokens[i][1], group, span
                        break

        for token, value, span in tokenize():
            # Handle partially parsed compound.
            if saved_token is not None and (
                    token in (_ReactionToken.Plus,
                              _ReactionToken.Arrow,
                              _ReactionToken.End)):
                compound = parse_compound(saved_token, global_comp)
                current_side.append((compound, 1))
                expect_operator = True
                saved_token = None

            if token == _ReactionToken.Plus:
                # Compound separator. Expect compound name or compound count
                # next.
                if not expect_operator:
                    raise ParseError('Unexpected token: {!r}'.format(value),
                                     span=span)
                expect_operator = False
            elif token == _ReactionToken.Arrow:
                # Reaction arrow. Expect compound name or compound count next.
                if direction is not None:
                    raise ParseError(
                        'More than one equation arrow: {!r}'.format(value),
                        span=span)

                if not expect_operator and len(left) > 0:
                    raise ParseError('Unexpected token: {!r}'.format(value),
                                     span=span)

                expect_operator = False
                direction = self._arrows[value]
                current_side = right
            elif token == _ReactionToken.Group:
                # Compound count. Expect compound name next.
                if expect_operator:
                    raise ParseError('Expected plus or arrow: {!r}'.format(
                        value), span=span)
                if saved_token is not None:
                    raise ParseError('Expected compound name: {!r}'.format(
                        value), span=span)
                saved_token = value
            elif token == _ReactionToken.Quoted:
                # Compound name. Expect operator next.
                if expect_operator:
                    raise ParseError('Expected plus or arrow: {!r}'.format(
                        value), span=span)
                if saved_token is not None:
                    try:
                        count = parse_compound_count(saved_token)
                    except ValueError as e:
                        raise_from(ParseError(
                            'Unable to parse compound count: {!r}'.format(
                                saved_token),
                            span=span), e)
                else:
                    count = 1
                compound = parse_compound(value, global_comp)
                current_side.append((compound, count))
                expect_operator = True
                saved_token = None
            elif token == _ReactionToken.Other:
                # Could be count or compound name. Store and expect other,
                # quoted, operator or end.
                if expect_operator:
                    raise ParseError('Expected plus or arrow: {!r}'.format(
                        value), span=span)
                if saved_token is not None:
                    try:
                        count = parse_compound_count(saved_token)
                    except ValueError as e:
                        raise_from(ParseError(
                            'Unable to parse compound count: {!r}'.format(
                                saved_token),
                            span=span), e)
                    compound = parse_compound(value, global_comp)
                    current_side.append((compound, count))
                    expect_operator = True
                    saved_token = None
                else:
                    saved_token = value

        if direction is None:
            raise ParseError('Missing equation arrow')

        return Reaction(direction, left, right)


_DEFAULT_PARSER = ReactionParser()


def parse_reaction(s):
    """Parse reaction string using the default parser."""
    return _DEFAULT_PARSER.parse(s)


def parse_compound(s, global_compartment=None):
    """Parse a compound specification.

    If no compartment is specified in the string, the global compartment
    will be used.
    """
    m = re.match(r'^\|(.*)\|$', s)
    if m:
        s = m.group(1)

    m = re.match(r'^(.+)\[(\S+)\]$', s)
    if m:
        compound_id = m.group(1)
        compartment = m.group(2)
    else:
        compound_id = s
        compartment = global_compartment

    return Compound(compound_id, compartment=compartment)


def parse_compound_count(s):
    """Parse a compound count (number of compounds)."""
    m = re.match(r'^\((.*)\)$', s)
    if m:
        s = m.group(1)

    for count_type in (int, Decimal, affine.Expression):
        try:
            return count_type(s)
        except:
            pass

    raise ValueError('Unable to parse compound count: {}'.format(s))
