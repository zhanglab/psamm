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
# Copyright 2019-2020       Jing Wang <jingwang89@uri.edu>

"""Importer for the COBRApy MAT format."""

from psamm.importer import Importer as BaseImporter, ModelLoadError
from psamm.datasource.entry import (DictCompoundEntry as CompoundEntry,
                                    DictReactionEntry as ReactionEntry,
                                    DictCompartmentEntry as CompartmentEntry)
from psamm.datasource import native
from psamm.reaction import Reaction, Compound, Direction
import re
import numpy as np
from scipy.io import loadmat
import decimal
import logging
import glob
import os


logger = logging.getLogger(__name__)


def _float_parser(num_str):
    num = float(num_str)
    if num.is_integer():
        return int(num)
    else:
        return decimal.Decimal(num_str)


class Importer(BaseImporter):
    """Read metabolic model from COBRA MATLAB mat format."""

    name = 'matlab'
    title = 'COBRA MATLAB mat'
    generic = True

    def help(self):
        """Print help text for importer."""
        print('Source must contain the model definition in COBRA MATLAB mat'
              ' format.\n'
              'Expected files in source directory:\n'
              '- *.mat')

    def _resolve_source(self, source):
        """Resolve source to filepath if it is a directory."""
        if os.path.isdir(source):
            sources = glob.glob(os.path.join(source, '*.mat'))
            if len(sources) == 0:
                raise ModelLoadError('No .mat file found in source directory')
            elif len(sources) > 1:
                raise ModelLoadError(
                    'More than one .mat file found in source directory')
            return sources[0]
        return source

    def _import(self, file, name):
        model_doc = loadmat(file)
        if name is None:
            names = [k for k in model_doc.keys()
                     if k not in ['__header__', '__version__', '__globals__']]
            modelnames = list()
            for i in names:
                if (isinstance(model_doc[i], np.ndarray) and
                        len(model_doc[i].dtype) > 1):
                    modelnames.append(i)
            if len(modelnames) == 1:
                name = str(modelnames[0])
                model_doc = model_doc[modelnames[0]]
            elif len(modelnames) > 1:
                raise ModelLoadError(
                    ('More than one model exist, '
                     'please choose from: %s '
                     'using --model-name' % modelnames))
            else:
                raise ModelLoadError('Incorrect format')
        else:
            if name not in model_doc.keys():
                raise ModelLoadError('Incorrect model name specified')
            else:
                model_doc = model_doc[name]
                if not (isinstance(model_doc, np.ndarray) and
                        len(model_doc.dtype) > 1):
                    raise ModelLoadError('Incorrect format')
        model = native.NativeModel()
        model.name = name
        self._compartments_from_compound = set()
        self._origid_compound = dict()
        model.compounds.update(self._read_compounds(model_doc))
        model.compartments.update(self._compartments_from_compound)

        self._boundary_compartment = None
        # Add model level compartment information
        if all(var in model_doc.dtype.names for var in ['comps', 'compNames']):
            model.compartments.update(self._read_compartments(model_doc))
            for comp in model.compartments:
                if comp.name.lower() == 'boundary':  # has boundary compartment
                    self._boundary_compartment = comp.id
            model.compartments.discard(self._boundary_compartment)

        for reaction, lower, upper in self._read_reactions(model_doc):
            model.reactions.add_entry(reaction)
            if lower is not None and upper is not None:
                model.limits[reaction.id] = reaction.id, lower, upper

        biomass_reaction = None
        objective_reactions = set()
        if 'c' in model_doc.dtype.names:
            for i in range(len(model_doc['c'][0, 0])):
                if model_doc['c'][0, 0][i][0] != 0:
                    objective_reactions.add(
                        str(model_doc['rxns'][0, 0][i][0][0]))
            if len(objective_reactions) == 1:
                biomass_reaction = next(iter(objective_reactions))
                logger.info('Detected biomass reaction: {}'.format(
                    biomass_reaction))
            elif len(objective_reactions) > 1:
                logger.warning(
                    'Multiple reactions are used as the'
                    ' biomass reaction: {}'.format(objective_reactions))
        else:
            logger.warning('No objective reaction')

        model.biomass_reaction = biomass_reaction

        return model

    def _read_compartments(self, doc):
        for i in range(len(doc['comps'][0, 0])):
            compartment = str(doc['comps'][0, 0][i][0][0])
            name = str(doc['compNames'][0, 0][i][0][0])
            yield CompartmentEntry(dict(id=compartment, name=name))

    def _read_compounds(self, doc):
        for i in range(len(doc['mets'][0, 0])):
            properties = dict()
            id = doc['mets'][0, 0][i][0][0]
            match = re.search(r'\[[0-9a-zA-Z]+\]$', id)  # e.g. h2o[c]
            uniqid = id
            compartment = None
            if match is not None:  # compartment information in id
                compartment = match.group().strip('[]')
                self._compartments_from_compound.add(
                    CompartmentEntry(dict(id=compartment, name='')))
                uniqid = re.sub(r'\[[0-9a-zA-Z]+\]$', '', id)
            else:
                match = re.search(r'_[0-9a-zA-Z]+$', id)  # e.g. h2o_c
                if match is not None:  # compartment information in id
                    compartment = match.group().strip('_')
                    self._compartments_from_compound.add(
                        CompartmentEntry(
                            dict(id=compartment, name='')))
                    uniqid = re.sub(r'_[0-9a-zA-Z]+$', '', id)
                else:
                    match = re.search(r'\([0-9a-zA-Z]+\)$', id)  # e.g. h2o(c)
                    if match is not None:  # compartment information in id
                        compartment = match.group().strip('()')
                        self._compartments_from_compound.add(
                            CompartmentEntry(dict(id=compartment, name='')))
                        uniqid = re.sub(r'\([0-9a-zA-Z]+\)$', '', id)
                    else:
                        match = re.search(r'[0-9][a-z]$', id)  # e.g. cpd100c
                        if match is not None:
                            compartment = str(match.group()[1])
                            self._compartments_from_compound.add(
                                CompartmentEntry(
                                    dict(id=compartment, name='')))
                            uniqid = id[:len(id) - 1]
            properties['id'] = uniqid
            self._origid_compound[id] = Compound(uniqid, compartment)
            for dtype in doc.dtype.names:
                if dtype.startswith('met') and dtype != 'mets':
                    # Compound related dtypes
                    if dtype == 'metFormulas':
                        if len(doc['metFormulas'][0, 0][i][0]) > 0:
                            properties[
                                'formula'] = doc['metFormulas'][0, 0][i][0][0]
                    elif dtype == 'metCharges':
                        properties['charge'] = _float_parser(
                            doc['metCharges'][0, 0][i][0])
                    elif dtype == 'metNames':
                        properties['name'] = str(
                            doc['metNames'][0, 0][i][0][0])
                    else:  # other dtypes
                        value = doc[dtype][0, 0][i][0]
                        if isinstance(value, np.ndarray):
                            if len(value) == 1:
                                value = value[0]
                            elif len(value) == 0:
                                continue
                        dtype = dtype.lstrip('met')
                        properties[dtype] = str(value)

            entry = CompoundEntry(properties)
            if 'formula' in entry.properties:
                formula = self._try_parse_formula(
                    entry.id, entry.properties['formula'])
                if formula is not None:
                    entry.properties['formula'] = formula
            yield entry

    def _read_reactions(self, doc):
        for i in range(len(doc['rxns'][0, 0])):
            properties = dict()
            properties['id'] = str(doc['rxns'][0, 0][i][0][0])

            for dtype in doc.dtype.names:
                if dtype.startswith('rxn') and dtype != 'rxns':
                    # reaction related dtypes
                    if dtype == 'rxnNames':
                        if len(doc['rxnNames'][0, 0][i][0]) > 0:
                            properties['name'] = str(
                                doc['rxnNames'][0, 0][i][0][0])
                    else:  # other dtypes
                        value = doc[dtype][0, 0][i][0]
                        if isinstance(value, np.ndarray):
                            if len(value) == 1:
                                value = value[0]
                            elif len(value) == 0:
                                continue
                        dtype = dtype.lstrip('rxn')
                        properties[dtype] = str(value)

            lower = float(doc['lb'][0, 0][i][0])
            upper = float(doc['ub'][0, 0][i][0])
            if lower == 0 and upper != 0:
                direction = Direction.Forward
            elif lower != 0 and upper == 0:
                direction = Direction.Reverse
            else:
                direction = Direction.Both

            # parse equation
            compounds = list()
            try:
                # incase it's a sparse matrix
                rows = doc['S'][0, 0][:, i].indices
            except AttributeError:
                # incase it's an ndarray
                rows = np.flatnonzero(doc['S'][0, 0][:, i])
            for j in rows:
                cpd = self._origid_compound.get(
                    doc['mets'][0, 0][j][0][0]
                )
                if cpd.compartment == self._boundary_compartment:
                    # skip compound that is out of boundary
                    continue
                compounds.append(
                    (
                        cpd,  # compound id
                        float(doc['S'][0, 0][j, i])  # stochiometry
                    )
                )
            properties['equation'] = (Reaction(direction, compounds))

            # parse genes
            if 'grRules' in doc.dtype.names:
                if len(doc['grRules'][0, 0][i][0]) > 0:
                    genes = doc['grRules'][0, 0][i][0][0]
                    if isinstance(genes, np.unicode_):
                        properties['genes'] = self._try_parse_gene_association(
                            properties['id'], genes
                        )

            # parse subsystem
            if 'subSystems' in doc.dtype.names:
                if doc['subSystems'][0, 0][i][0].size > 0:
                    # the subSystem has either 2 levels or 4 levels
                    # check which level to obtain
                    if isinstance(doc['subSystems'][0, 0][i][0][0], str):
                        properties['subsystem'] = str(
                            doc['subSystems'][0, 0][i][0][0])
                    elif len(doc['subSystems'][0, 0][i][0][0][0]) > 0:
                        properties['subsystem'] = str(
                            doc['subSystems'][0, 0][i][0][0][0][0])

            entry = ReactionEntry(properties)
            yield entry, lower, upper

    def import_model(self, source, model_name=None):
        """Import and return model instance."""
        if not hasattr(source, 'read'):  # Not a File-like object
            with open(self._resolve_source(source), 'rb') as f:
                return self._import(f, model_name)
        else:
            return self._import(source, model_name)
