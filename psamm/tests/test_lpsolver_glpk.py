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
# Copyright 2026  Jason Vailionis <jason_vailionis@uri.edu>

import resource
import subprocess
import sys
import textwrap
import unittest

try:
    from psamm.lpsolver import glpk
except ImportError:
    glpk = None

requires_solver = unittest.skipIf(glpk is None, 'GLPK solver not available')

# Solve a problem on the main thread and let the cyclic garbage collector
# free it on another thread. GLPK aborts the process if the problem is
# deleted from a thread other than the one that created it.
_CROSS_THREAD_DELETE = textwrap.dedent('''
    import gc
    import threading
    from psamm.lpsolver import glpk

    gc.disable()
    prob = glpk.Solver().create_problem()
    prob.define('x', lower=0, upper=10)
    prob.set_objective(prob.var('x'))
    prob.solve()
    del prob

    thread = threading.Thread(target=gc.collect)
    thread.start()
    thread.join()
    print('pending', sum(len(p) for p in glpk._pending_deletes.values()))

    # The deferred problem is deleted when the main thread uses GLPK again
    glpk.Solver().create_problem()
    print('pending', sum(len(p) for p in glpk._pending_deletes.values()))
''')


def _disable_core_dump():
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))


@requires_solver
class TestGLPKProblemDelete(unittest.TestCase):
    def test_delete_from_other_thread(self):
        # Run in a subprocess since a failure aborts the whole process
        proc = subprocess.run(
            [sys.executable, '-c', _CROSS_THREAD_DELETE],
            capture_output=True, text=True, timeout=60,
            preexec_fn=_disable_core_dump)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(proc.stdout.split('\n'), [
            'pending 1', 'pending 0', ''])


if __name__ == '__main__':
    unittest.main()
