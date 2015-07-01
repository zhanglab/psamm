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

"""Miscellaneaus data source functions."""

import re
from decimal import Decimal

from psamm.reaction import Reaction, Compound


def parse_sudensimple_reaction(s, arrow_rev='<=>', arrow_irrev='->'):
    """Parse a reaction string (SudenSimple)"""

    # Compile pattern for matching reaction arrows
    arrow_pattern = re.compile(
        '(' + ('|'.join(re.escape(x) for x in (arrow_irrev, arrow_rev))) + ')')

    def parse_compound_list(s):
        if s.strip() == '':
            return

        for cpd in s.strip().split('+'):
            cpd_split = cpd.strip().split(' ')
            if len(cpd_split) == 1:
                count = 1
                spec = cpd_split[0].strip()
            else:
                count = Decimal(cpd_split[0])
                if count % 1 == 0:
                    count = int(count)
                spec = cpd_split[1].strip()

            m = re.match(r'^(.+)\[(.+)\]$', spec)
            if m:
                # compartment
                cpdid = m.group(1)
                comp = m.group(2)
            else:
                cpdid = spec
                comp = None

            yield Compound(cpdid, compartment=comp), count

    # Split by equation arrow
    direction = Reaction.Right
    left, arrow, right = arrow_pattern.split(s, maxsplit=1)
    direction = Reaction.Right if arrow == arrow_irrev else Reaction.Bidir

    return Reaction(direction, parse_compound_list(left),
                    parse_compound_list(right))


def parse_metnet_reaction(s, arrow_rev='<==>', arrow_irrev='-->'):
    """Parser for the reaction format in MetNet model"""

    # Compile pattern for matching reaction arrows
    arrow_pattern = re.compile(
        '(' + ('|'.join(re.escape(x) for x in (arrow_irrev, arrow_rev))) + ')')

    def parse_compound_list(s, global_comp):
        if s.strip() == '':
            return

        for cpd in s.strip().split('+'):
            cpd_split = cpd.strip().split(') ')
            if len(cpd_split) == 1:
                count = 1
                spec = cpd_split[0].strip()
            else:
                count = Decimal(cpd_split[0].lstrip('('))
                if count % 1 == 0:
                    count = int(count)
                spec = cpd_split[1].strip()

            if global_comp is None:
                m = re.match(r'^(.+)\[(.+)\]$', spec)
                if m:
                    # compartment
                    cpdid = m.group(1)
                    comp = m.group(2)
                else:
                    cpdid = spec
                    comp = None
            else:
                cpdid = spec
                comp = global_comp

            yield Compound(cpdid, compartment=comp), count

    # Split by colon for compartment information
    m = re.match(r'^\s*\[(\w+)\]\s*:\s*(.*)', s)
    if m is not None:
        global_comp = m.group(1)
        s = m.group(2)
    else:
        global_comp = None

    # Split by equation arrow
    direction = Reaction.Right
    left, arrow, right = arrow_pattern.split(s, maxsplit=1)
    direction = Reaction.Right if arrow == arrow_irrev else Reaction.Bidir

    return Reaction(direction, parse_compound_list(left, global_comp),
                    parse_compound_list(right, global_comp))
