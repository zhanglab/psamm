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
# Copyright 2014-2017  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2018-2018  Ke Zhang <kzhang@my.uri.edu>

from __future__ import unicode_literals

import time
import logging
from ..command import (LoopRemovalMixin, ObjectiveMixin, SolverCommandMixin,
                       MetabolicMixin, Command)
import csv
from ..reaction import Direction
form six import text_type, iteritems
from .. import findprimatypairs
from ..formula import Formula, Atom, ParseError
from .. import graph
from collections import Counter


logger = logging.getLogger(__name__)

REACTION_COLOR = '#ccebc5'
COMPOUND_COLOR = '#b3cde3'
ACTIVE_COLOR = '#fbb4ae'
ALT_COLOR = '#ccb460'


