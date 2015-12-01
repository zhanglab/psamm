#!/usr/bin/env python
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
# Copyright 2015  Jon Lund Steffensen <jon_steffensen@uri.edu>

import os
import unittest

from psamm.lpsolver import generic


class TestListSolversCommand(unittest.TestCase):
    def test_list_lpsolvers(self):
        if os.environ.get('PSAMM_SOLVER') in ('nosolver', None):
            with self.assertRaises(SystemExit):
                generic.list_solvers([])
        else:
            generic.list_solvers([])
