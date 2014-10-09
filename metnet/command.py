
'''Utilities for the command line interface'''

import argparse

from metnet.metabolicmodel import DictDatabase

class Command(object):
    '''Represents a command in the interface, operating on a model or database

    Subclasses must implement __call__ to handle command execution. Arguments
    will be given as keyword arguments. The keywords database, model and
    compounds will be set.

    In addition, init_parser() can be implemented which will allow the
    command to initialize an instance of ArgumentParser as desired. The
    resulting arguments will be given as keyword arguments to __call__.
    '''

    def __init__(self, title=None):
        self._title = title

    @property
    def title(self):
        '''Command title'''
        return self._title

    def init_parser(self, parser):
        '''Initialize command line parser (argparse.ArgumentParser)'''
        pass

    def __call__(self):
        '''Execute command'''
        pass


def main(command):
    '''Run the command line interface with the given Command'''

    parser = argparse.ArgumentParser(description=command.title)
    parser.add_argument('--database', metavar='file', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Files to use as reaction database')
    parser.add_argument('--compounds', metavar='file', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Files to use as compound database')
    parser.add_argument('--model', metavar='model', nargs=1,
                        type=argparse.FileType('r'),
                        help='File to use as model definition (database subset)')

    command.init_parser(parser)
    args = parser.parse_args()

    # Load reaction database from file
    database = DictDatabase.load_from_files(*args.database)
    model = None
    if args.model is not None:
        # Set database and model to the database subset
        database = database.load_model_from_file(args.model[0])
        model = database

    compounds = args.compounds

    # Call command
    kwargs = { key: value for key, value in vars(args).iteritems() if key not in ('database', 'compounds', 'model') }
    command(database=database, model=model, compounds=compounds, **kwargs)