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
# Copyright 2018-2020  Jing Wang <wjingsjtu@gmail.com>
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>
# Copyright 2020-2021  Elysha Sameth <esameth1@uri.edu>

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
import errno
import pandas as pd
import re
from future import standard_library
standard_library.install_aliases()


class Curator(object):
    """Parse and save mapping files during manual curation.

    Use :meth:`.add_mapping` to add new curated pairs. Save current progress
    into files by :meth:`.save`.

    Besides the curated mapping files, the :class:`Curator` will also store
    false mappings into `.false` files, and compounds and reactions
    to be ignored can be stored
    in `.ignore` files. For example, if the `curated_compound_map_file` is set
    to `curated_compound_mapping.tsv`, then the false mappings will be stored
    in `curated_compound_mapping.tsv.false`, and the pairs to be ignored
    should be stored in `curated_compound_mapping.tsv.ignore`.

    If the curated files already exist, the :class:`Curator` will consider them
    as the previous progress, then append new curation results.

    Args:
        compound_map_file: .tsv file of compound mapping result
        reaction_map_file: .tsv file of reaction mapping result
        curated_compound_map_file: .tsv file of curated compound mapping result
        curated_reaction_map_file: .tsv file of curated reaction mapping result
    """

    def __init__(self,
                 compound_map_file, reaction_map_file,
                 curated_compound_map_file, curated_reaction_map_file):
        self.compound_map_file = compound_map_file
        self.reaction_map_file = reaction_map_file
        self.curated_compound_map_file = curated_compound_map_file
        self.curated_reaction_map_file = curated_reaction_map_file
        self.false_compound_map_file = curated_compound_map_file + '.false'
        self.false_reaction_map_file = curated_reaction_map_file + '.false'
        self.ignore_compound_file = curated_compound_map_file + '.ignore'
        self.ignore_reaction_file = curated_reaction_map_file + '.ignore'

        self.compound_map = read_mapping(self.compound_map_file, [0, 1])
        self.compound_map.sort_values(by='p', inplace=True, ascending=False)
        self.reaction_map = read_mapping(self.reaction_map_file, [0, 1])
        self.reaction_map.sort_values(by='p', inplace=True, ascending=False)
        self.curated_compound_map = read_mapping(
            self.curated_compound_map_file, [0, 1])
        self.curated_reaction_map = read_mapping(
            self.curated_reaction_map_file, [0, 1])
        self.false_compound_map = read_mapping(
            self.false_compound_map_file, [0, 1])
        self.false_reaction_map = read_mapping(
            self.false_reaction_map_file, [0, 1])
        self.ignore_compound = read_ignore(self.ignore_compound_file)
        self.ignore_reaction = read_ignore(self.ignore_reaction_file)

        self.num_curated_compounds = len(self.curated_compound_map.index)
        self.num_curated_compounds_left = \
            self.compound_map.index.get_level_values(0).nunique() - \
            self.num_curated_compounds
        self.num_curated_reactions = len(self.curated_reaction_map.index)
        self.num_curated_reactions_left = \
            self.reaction_map.index.get_level_values(0).nunique() - \
            self.num_curated_reactions

    def reaction_checked(self, id):
        """Return True if reaction pair has been checked.

        Args:
            id: one reaction id or a tuple of id pair
        """
        return (id in self.curated_reaction_map.index or
                id in self.false_reaction_map.index)

    def reaction_ignored(self, id):
        """Return True if reaction id is in ignore list."""
        return id in self.ignore_reaction

    def compound_checked(self, id):
        """Return True if compound pair has been checked.

        Args:
            id: one compound id or a tuple of id pair
        """
        return (id in self.curated_compound_map.index or
                id in self.false_compound_map.index)

    def compound_ignored(self, id):
        """Return True if compound id is in ignore list."""
        return id in self.ignore_compound

    def add_mapping(self, id, type, correct):
        """Add new mapping result to curated list.

        Args:
            id: tuple of id pair
            type: 'c' for compound mapping, 'r' for reaction mapping
            correct: True if the mapping pair is correct
        """
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

    def add_ignore(self, id, type):
        """Add id to ignore list.

        Args:
            id: id to be ignored
            type: 'c' for compound mapping, 'r' for reaction mapping
        """
        if type == 'c':
            self.ignore_compound.append(id)
        if type == 'r':
            self.ignore_reaction.append(id)

    def save(self):
        """Save current curator to files."""
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
        if len(self.ignore_compound) > 0:
            write_ignore(self.ignore_compound, self.ignore_compound_file)
        if len(self.ignore_reaction) > 0:
            write_ignore(self.ignore_reaction, self.ignore_reaction_file)
        print('Progress saved\n')


def read_mapping(file, index_col):
    try:
        df = pd.read_csv(file, sep='\t', index_col=index_col)
    except IOError as e:
        if e.errno != errno.ENOENT:
            raise
        df = pd.DataFrame()
    return df


def add_mapping(map, row):
    return map.append(pd.DataFrame(row).T)


def read_ignore(file):
    ignore = list()
    try:
        with open(file) as f:
            for r in f:
                ignore.append(r.strip())
    except IOError as e:
        if e.errno != errno.ENOENT:
            raise
    return ignore


def write_ignore(ignore, file):
    with open(file, 'w') as o:
        for i in ignore:
            o.write('%s\n' % str(i))


def filter_search_term(s):
    return re.sub(r'[^a-z0-9]+', '', s.lower())


def search_compound(model, id):
    """Search a set of compounds, then print detailed properties.

    Args:
        id: a list of compound ids
    """
    selected_compounds = set()

    for compound in model.compounds:
        if len(id) > 0:
            if any(c == compound.id for c in id):
                selected_compounds.add(compound)
                continue

    # Show results
    for compound in selected_compounds:
        props = set(compound.properties) - {'id'}
        if compound.filemark is not None:
            print('Defined in {}'.format(compound.filemark))
        print('id: {}'.format(compound.id))
        for prop in sorted(props):
            print('{}: {}'.format(prop, compound.properties[prop]))
        print('\n')


def search_reaction(model, ids):
    """Search a set of reactions, print detailed properties, then return a
    generator. Each item in the generator is a list of compounds in the
    corresponding reaction.

    Args:
        ids: a list of reaction ids
    """
    selected_reactions = set()
    for r in ids:
        if r in model.reactions:
            selected_reactions.add(model.reactions[r])

    # Show results
    for reaction in selected_reactions:
        props = set(reaction.properties) - {'id', 'equation'}
        if reaction.filemark is not None:
            print('Defined in {}'.format(reaction.filemark))
        print('id: {}'.format(reaction.id))
        print('equation: {}'.format(
            reaction.equation))
        translated_equation = reaction.equation.translated_compounds(
            lambda c: model.compounds[c].name)
        if reaction.equation != translated_equation:
            print('equation (compound names): {}'.format(
                translated_equation))
        for prop in sorted(props):
            print('{}: {}'.format(prop, reaction.properties[prop]))
        print('\n')
        yield [c for c, v in reaction.equation.compounds]
