#!/usr/bin/env python

import os
import logging

# Set up logging for the command line interface
# This must be set up before any of the module-level
# loggers are instantiated.
if 'DEBUG' in os.environ:
    level = getattr(logging, os.environ['DEBUG'].upper(), None)
    if level is not None:
        logging.basicConfig(level=level)
else:
    logging.basicConfig(level=logging.INFO)

from metnet import command

if __name__ == '__main__':
    command.main()
