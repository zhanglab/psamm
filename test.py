#!/usr/bin/env python

import unittest

if __name__ == '__main__':
    # Run unit tests in tests directory
    testsuite = unittest.TestLoader().discover('metnet/tests')
    unittest.TextTestRunner(verbosity=2).run(testsuite)
