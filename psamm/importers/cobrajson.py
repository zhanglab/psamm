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

"""Importer for the COBRApy JSON format."""

import os
import glob
import json
import logging
import decimal

from six import iteritems

from ..reaction import Reaction, Compound, Direction
from ..datasource import native
from ..datasource.entry import (DictCompoundEntry as CompoundEntry,
                                DictReactionEntry as ReactionEntry,
                                DictCompartmentEntry as CompartmentEntry)

from ..importer import Importer as BaseImporter, ModelLoadError

logger = logging.getLogger(__name__)


def _float_parser(num_str):
    num = float(num_str)
    if num.is_integer():
        return int(num)
    else:
        return decimal.Decimal(num_str)


class Importer(BaseImporter):
    """Read metabolic model from COBRApy JSON format."""

    name = 'json'
    title = 'COBRApy JSON'
    generic = True

    def help(self):
        """Print help text for importer."""
        print('Source must contain the model definition in COBRApy JSON'
              ' format.\n'
              'Expected files in source directory:\n'
              '- *.json')

    def _resolve_source(self, source):
        """Resolve source to filepath if it is a directory."""
        if os.path.isdir(source):
            sources = glob.glob(os.path.join(source, '*.json'))
            if len(sources) == 0:
                raise ModelLoadError('No .json file found in source directory')
            elif len(sources) > 1:
                raise ModelLoadError(
                    'More than one .json file found in source directory')
            return sources[0]
        return source

    def _import(self, file):
        model_doc = json.load(file, parse_float=_float_parser)
        model = native.NativeModel()
        model.compartments.update(self._read_compartments(model_doc))
        model.compounds.update(self._read_compounds(model_doc))

        for reaction, lower, upper in self._read_reactions(model_doc):
            model.reactions.add_entry(reaction)
            if lower is not None and upper is not None:
                model.limits[reaction.id] = reaction.id, lower, upper

        model.name = model_doc.get('id', 'COBRA JSON model')

        biomass_reaction = None
        objective_reactions = set()
        for reaction in model_doc['reactions']:
            if reaction.get('objective_coefficient', 0) != 0:
                objective_reactions.add(reaction['id'])

        if len(objective_reactions) == 1:
            biomass_reaction = next(iter(objective_reactions))
            logger.info('Detected biomass reaction: {}'.format(
                biomass_reaction))
        elif len(objective_reactions) > 1:
            logger.warning(
                'Multiple reactions are used as the'
                ' biomass reaction: {}'.format(objective_reactions))

        model.biomass_reaction = biomass_reaction

        # Set compartment in reaction compounds
        compound_compartment = {}
        for compound in model.compounds:
            compartment = compound.properties.pop('compartment')
            compound_compartment[compound.id] = compartment

        for reaction in model.reactions:
            equation = reaction.equation
            if equation is None:
                continue

            # Translate compartments in reaction
            left = ((c.in_compartment(compound_compartment[c.name]), v) for
                    c, v in equation.left)
            right = ((c.in_compartment(compound_compartment[c.name]), v) for
                     c, v in equation.right)
            reaction.properties['equation'] = Reaction(
                equation.direction, left, right)

        return model

    def _read_compartments(self, doc):
        for compartment, name in iteritems(doc.get('compartments', {})):
            yield CompartmentEntry(dict(id=compartment, name=name))

    def _read_compounds(self, doc):
        for compound in doc['metabolites']:
            entry = CompoundEntry(compound)
            if 'formula' in entry.properties:
                formula = self._try_parse_formula(
                    entry.id, entry.properties['formula'])
                if formula is not None:
                    entry.properties['formula'] = formula
            yield entry

    def _parse_reaction_equation(self, entry):
        metabolites = entry.properties.pop('metabolites')
        compounds = ((Compound(metabolite), value)
                     for metabolite, value in iteritems(metabolites))
        if (entry.properties.get('lower_bound') == 0 and
                entry.properties.get('upper_bound') != 0):
            direction = Direction.Forward
        elif (entry.properties.get('lower_bound') != 0 and
              entry.properties.get('upper_bound') == 0):
            direction = Direction.Reverse
        else:
            direction = Direction.Both
        return Reaction(direction, compounds)

    def _read_reactions(self, doc):
        for reaction in doc['reactions']:
            entry = ReactionEntry(reaction)

            entry.properties['equation'] = (
                self._parse_reaction_equation(entry))

            lower, upper = None, None

            if 'lower_bound' in entry.properties:
                lower = entry.properties.pop('lower_bound')
            if 'upper_bound' in entry.properties:
                upper = entry.properties.pop('upper_bound')

            if 'gene_reaction_rule' in entry.properties:
                genes = self._try_parse_gene_association(
                    entry.id, entry.properties.pop('gene_reaction_rule'))
                if genes is not None:
                    entry.properties['genes'] = genes

            yield entry, lower, upper

    def import_model(self, source):
        """Import and return model instance."""
        if not hasattr(source, 'read'):  # Not a File-like object
            with open(self._resolve_source(source), 'r') as f:
                return self._import(f)
        else:
            return self._import(source)
