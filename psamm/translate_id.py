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
# Copyright 2018-2019  Jing Wang <wjingsjtu@gmail.com>

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from builtins import open
from builtins import super
from future import standard_library
standard_library.install_aliases()  # noqa
import csv
import logging
from psamm.datasource.native import NativeModel
from psamm.importer import write_yaml_model


def read_mapping(mapping):
    mapping_file = open(mapping, mode='rU')
    mapping_id = {}
    for row in csv.reader(mapping_file, delimiter='\t'):
        if row[1] == '':
            continue
        else:
            mapping_id[row[0]] = row[1]
    return(mapping_id)


def translate_id(id, mapping):
    if id not in list(mapping):
        return id
    return mapping[id]


class TranslatedModel(NativeModel):
    """A NativeModel with translated ids based on reference model.

    The compound_map and reaction_map are two :class:`dict` that use original
    id as key and new id as value.
    Use :meth:`write_model` to output yaml files.
    """
    def __init__(self, compound_map, reaction_map, ref_model):
        super().__init__()
        self._cpd_mapping_id = compound_map
        self._rxn_mapping_id = reaction_map

        # initial a logger
        self._logger = logging.getLogger(__name__)

        # inherit properties
        self._properties = ref_model._properties
        self.version_string = 'model mapped from: %s' % self.version_string
        self._compartment_boundaries = ref_model.compartment_boundaries
        self._compartments = ref_model.compartments

        # translate ids
        self._translate_compounds(ref_model)
        self._translate_reactions(ref_model)
        self._translate_exchange(ref_model)
        self._translate_limits(ref_model)

    @property
    def cpd_mapping_id(self):
        return self._cpd_mapping_id

    @property
    def rxn_mapping_id(self):
        return self._rxn_mapping_id

    def _translate_compounds(self, ref_model):
        """Translate compound ids and corresponding exchange reactions."""
        for cpd in ref_model.compounds:
            new_id = translate_id(cpd.id, self.cpd_mapping_id)
            if new_id not in self.compounds:  # do not duplicate compounds
                cpd.properties['id'] = new_id
                cpd.properties['original_id'] = [cpd.id]
                # refresh with changed properties
                cpd = cpd.__class__(cpd.properties, filemark=cpd.filemark)
                self.compounds.add_entry(cpd)
            else:
                self.compounds[new_id].properties['original_id'].append(cpd.id)
                self._logger.warning(
                    'Multiple compounds: %s will be translated to %s',
                    self.compounds[new_id].properties['original_id'],
                    new_id)

    def _compound_trans(self, id):
        """Make a directly callable method for compound.translate()"""
        return translate_id(id, self.cpd_mapping_id)

    def _translate_reactions(self, ref_model):
        """Translate reaction ids, equations and limits."""

        for rxn in ref_model.reactions:
            new_id = translate_id(rxn.id, self.rxn_mapping_id)
            if new_id not in self.reactions:  # do not duplicate reactions
                rxn.properties['id'] = new_id
                rxn.properties['original_id'] = rxn.id
                rxn.equation = rxn.equation.translated_compounds(
                    self._compound_trans)
                # test whether the same compound occurs at both sides
                compound_common = set(rxn.equation.left).intersection(
                    set(rxn.equation.right))
                if len(compound_common) > 0:
                    self._logger.error((
                        'Reaction %s will have '
                        'the same compounds %s on both sides: %s'),
                        rxn.id, compound_common, rxn.equation)
                # refresh with changed properties
                rxn = rxn.__class__(rxn.properties, filemark=rxn.filemark)
                self.reactions.add_entry(rxn)
            else:
                self._logger.error(
                    ('Reaction %s is omitted because its new id %s has been '
                     'assigned to %s'),
                    rxn.id,
                    new_id,
                    self.reactions[new_id].properties['original_id'])

    def _translate_exchange(self, ref_model):
        """Translate exchange reactions."""
        # cpd is :class:`Compound(name, compartment)`
        for cpd in ref_model.exchange:
            new_cpd = cpd.translate(self._compound_trans)
            if new_cpd not in self.exchange:  # do not duplicate exchange
                lower, upper = ref_model.exchange[cpd][2:]
                name = 'EX_%s(e)' % new_cpd.name
                self.exchange[new_cpd] = (new_cpd, name, lower, upper)
            else:
                self._logger.error(
                    ('Exchange reaction of compound %s is omitted because '
                     'its new id %s has been assigned multiple compounds: %s'),
                    cpd.name,
                    new_cpd.name,
                    self.compounds[new_cpd.name].properties['original_id'])

    def _translate_limits(self, ref_model):
        """"Translate reaction limits."""
        for rxn in ref_model.limits:
            new_rxn = translate_id(rxn, self.rxn_mapping_id)
            if new_rxn not in self.limits:  # do not duplicate limits
                lower, upper = ref_model.limits[rxn][1:]
                self.limits[new_rxn] = (new_rxn, lower, upper)
            else:
                self._logger.error(
                    ('Limits of reaction %s is omitted because its new id %s '
                     'has been assigned to %s'),
                    rxn,
                    new_rxn,
                    self.reactions[new_rxn].properties['original_id'])

    def write_model(self, dest, **kwargs):
        """Output model into YAML files."""
        write_yaml_model(self, dest, **kwargs)
