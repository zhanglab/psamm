from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
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

import errno
import pandas as pd
import re
from future import standard_library
standard_library.install_aliases()


class Curator(object):

    def __init__(self,
                 compound_map_file, reaction_map_file,
                 curated_compound_map_file, curated_reaction_map_file):
        self.compound_map_file = compound_map_file
        self.reaction_map_file = reaction_map_file
        self.curated_compound_map_file = curated_compound_map_file
        self.curated_reaction_map_file = curated_reaction_map_file
        self.false_compound_map_file = curated_compound_map_file + '.false'
        self.false_reaction_map_file = curated_reaction_map_file + '.false'
        self.compound_map = read_mapping(self.compound_map_file)
        self.compound_map.sort_values(by='p', inplace=True, ascending=False)
        self.reaction_map = read_mapping(self.reaction_map_file)
        self.reaction_map.sort_values(by='p', inplace=True, ascending=False)
        self.curated_compound_map = read_mapping(
            self.curated_compound_map_file)
        self.curated_reaction_map = read_mapping(
            self.curated_reaction_map_file)
        self.false_compound_map = read_mapping(
            self.false_compound_map_file)
        self.false_reaction_map = read_mapping(
            self.false_reaction_map_file)

    def reaction_checked(self, id):
        return (id in self.curated_reaction_map.index or
                id in self.false_reaction_map.index)

    def compound_checked(self, id):
        return (id in self.curated_compound_map.index or
                id in self.false_compound_map.index)

    def add_mapping(self, id, type, correct):
        if type == 'c':
            if correct:
                self.curated_compound_map = add_mapping(
                    self.curated_compound_map,
                    self.compound_map.loc[id]
                )
            else:
                self.false_compound_map = add_mapping(
                    self.false_compound_map,
                    self.compound_map.loc[id]
                )
        if type == 'r':
            if correct:
                self.curated_reaction_map = add_mapping(
                    self.curated_reaction_map,
                    self.reaction_map.loc[id]
                )
            else:
                self.false_reaction_map = add_mapping(
                    self.false_reaction_map,
                    self.reaction_map.loc[id]
                )

    def save(self):
        if len(self.curated_compound_map) > 0:
            self.curated_compound_map.to_csv(
                self.curated_compound_map_file, sep='\t'
            )
        if len(self.curated_reaction_map) > 0:
            self.curated_reaction_map.to_csv(
                self.curated_reaction_map_file, sep='\t'
            )
        if len(self.false_compound_map) > 0:
            self.false_compound_map.to_csv(
                self.false_compound_map_file, sep='\t'
            )
        if len(self.false_reaction_map) > 0:
            self.false_reaction_map.to_csv(
                self.false_reaction_map_file, sep='\t'
            )
        print('Progress saved\n')


def read_mapping(file):
    try:
        df = pd.read_csv(file, sep='\t', index_col=[0, 1])
    except OSError as e:
        if e.errno != errno.ENOENT:
            raise
        df = pd.DataFrame()
    return df


def add_mapping(map, row):
    return map.append(pd.DataFrame(row).T)


def filter_search_term(s):
    return re.sub(r'[^a-z0-9]+', '', s.lower())


def search_compound(model, id):
    selected_compounds = set()

    for compound in model.compounds:
        if len(id) > 0:
            if any(c == compound.id for c in id):
                selected_compounds.add(compound)
                continue

    # Show results
    for compound in selected_compounds:
        props = set(compound.properties) - {'id'}
        print('id: {}'.format(compound.id))
        for prop in sorted(props):
            print('{}: {}'.format(prop, compound.properties[prop]))
        if compound.filemark is not None:
            print('Defined in {}'.format(compound.filemark))
        print('\n')


def search_reaction(model, id):
    selected_reactions = set()

    # Prepare translation table from compound id to name
    compound_name = {}
    for compound in model.compounds:
        if 'name' in compound.properties:
            compound_name[compound.id] = compound.properties['name']

    for reaction in model.reactions:
        if len(id) > 0:
            if any(r == reaction.id for r in id):
                selected_reactions.add(reaction)
                continue

    # Show results
    for reaction in selected_reactions:
        props = set(reaction.properties) - {'id', 'equation'}
        print('id: {}'.format(reaction.id))
        print('equation: {}'.format(
            reaction.equation))
        translated_equation = reaction.equation.translated_compounds(
            lambda x: compound_name.get(x, x))
        if reaction.equation != translated_equation:
            print('equation (compound names): {}'.format(
                translated_equation))
        for prop in sorted(props):
            print('{}: {}'.format(prop, reaction.properties[prop]))
        if reaction.filemark is not None:
            print('Defined in {}'.format(reaction.filemark))
        print('\n')
        return [c.name for c, v in reaction.equation.compounds]
