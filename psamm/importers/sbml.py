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
# Copyright 2015-2017  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>

"""SBML importers."""

import os
import glob
import logging


from ..datasource import sbml
from ..datasource.context import FilePathContext
from ..datasource.entry import (DictCompartmentEntry, DictCompoundEntry,
                                DictReactionEntry)

from ..importer import Importer, ParseError, ModelLoadError

logger = logging.getLogger(__name__)


class BaseImporter(Importer):
    """Base importer for reading metabolic model from an SBML file."""

    def help(self):
        """Print importer help text."""
        print('Source must contain the model definition in SBML format.\n'
              'Expected files in source directory:\n'
              '- *.sbml')

    def _resolve_source(self, source):
        """Resolve source to filepath if it is a directory."""
        if os.path.isdir(source):
            sources = glob.glob(os.path.join(source, '*.sbml'))
            if len(sources) == 0:
                raise ModelLoadError('No .sbml file found in source directory')
            elif len(sources) > 1:
                raise ModelLoadError(
                    'More than one .sbml file found in source directory')
            return sources[0]
        return source

    def _get_reader(self, f):
        raise NotImplementedError('Subclasses must implement _get_reader()')

    def import_model(self, source):
        """Import and return model instance."""
        source = self._resolve_source(source)
        self._context = FilePathContext(source)
        with self._context.open() as f:
            self._reader = self._open_reader(f)

        return self._reader.create_model()


class StrictImporter(BaseImporter):
    """Read metabolic model from an SBML file using strict parser."""

    name = 'SBML-strict'
    title = 'SBML model (strict)'
    generic = True

    def _open_reader(self, f):
        try:
            return sbml.SBMLReader(f, strict=True, ignore_boundary=True)
        except sbml.ParseError as e:
            raise ParseError(e)

    def _translate_compartment(self, entry, new_id):
        return DictCompartmentEntry(entry, new_id)

    def _translate_compound(self, entry, new_id, compartment_map):
        return DictCompoundEntry(entry, new_id)

    def _translate_reaction(
            self, entry, new_id, compartment_map, compound_map):
        return DictReactionEntry(entry, new_id)

    def import_model(self, source):
        """Import and return model instance."""
        model = super(StrictImporter, self).import_model(source)

        # Translate entries into dict-based entries.
        sbml.convert_model_entries(
            model, convert_id=lambda entry: entry.id,
            translate_compartment=self._translate_compartment,
            translate_compound=self._translate_compound,
            translate_reaction=self._translate_reaction)

        return model


class NonstrictImporter(BaseImporter):
    """Read metabolic model from an SBML file using non-strict parser."""

    name = 'SBML'
    title = 'SBML model (non-strict)'
    generic = True

    def _open_reader(self, f):
        try:
            return sbml.SBMLReader(f, strict=False, ignore_boundary=True)
        except sbml.ParseError as e:
            raise ParseError(e)

    def import_model(self, source):
        """Import and return model instance."""
        model = super(NonstrictImporter, self).import_model(source)
        sbml.convert_sbml_model(model)
        return model
