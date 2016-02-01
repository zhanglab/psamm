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

from .reaction import ReactionParser
from ..reaction import Direction


def parse_sudensimple_reaction(s, arrow_rev='<=>', arrow_irrev='->'):
    """Parse a reaction string (SudenSimple).

    .. deprecated:: 0.18
       Use :func:`~psamm.datasource.reaction.parse_reaction` instead.
    """

    arrows = (
        (arrow_rev, Direction.Both),
        (arrow_irrev, Direction.Forward)
    )

    return ReactionParser(arrows=arrows).parse(s)


def parse_metnet_reaction(s, arrow_rev='<==>', arrow_irrev='-->'):
    """Parser for the reaction format in MetNet model.

    .. deprecated:: 0.18
       Use :func:`~psamm.datasource.reaction.parse_reaction` instead.
    """

    arrows = (
        (arrow_irrev, Direction.Forward),
        (arrow_rev, Direction.Both)
    )

    return ReactionParser(arrows=arrows, parse_global=True).parse(s)
