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
# Copyright 2014-2017  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2016  Chao Liu <lcddzyx@gmail.com>
# Copyright 2020 Christopher Powers <c-11060@my.uri.edu>
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>

"""Functionality related to gap-filling in general.

This module contains some general functions for preparing models for
gap-filling. Specific gap-filling methods are implemented in the ``gapfill``
and ``fastgapfill`` modules.
"""

from __future__ import unicode_literals

import logging

from six import iteritems

from .metabolicmodel import create_exchange_id, create_transport_id
from .reaction import Reaction, Direction

logger = logging.getLogger(__name__)


def add_all_database_reactions(model, compartments):
    """Add all reactions from database that occur in given compartments.

    Args:
        model: :class:`psamm.metabolicmodel.MetabolicModel`.
    """

    added = set()
    for rxnid in model.database.reactions:
        reaction = model.database.get_reaction(rxnid)
        if all(compound.compartment in compartments
               for compound, _ in reaction.compounds):
            if not model.has_reaction(rxnid):
                added.add(rxnid)
            model.add_reaction(rxnid)

    return added


def add_all_exchange_reactions(model, compartment, allow_duplicates=False):
    """Add all exchange reactions to database and to model.

    Args:
        model: :class:`psamm.metabolicmodel.MetabolicModel`.
    """

    all_reactions = {}
    if not allow_duplicates:
        # TODO: Avoid adding reactions that already exist in the database.
        # This should be integrated in the database.
        for rxnid in model.database.reactions:
            rx = model.database.get_reaction(rxnid)
            all_reactions[rx] = rxnid

    added = set()
    added_compounds = set()
    initial_compounds = set(model.compounds)
    reactions = set(model.database.reactions)
    for model_compound in initial_compounds:
        compound = model_compound.in_compartment(compartment)
        if compound in added_compounds:
            continue

        rxnid_ex = create_exchange_id(reactions, compound)

        reaction_ex = Reaction(Direction.Both, {compound: -1})
        if reaction_ex not in all_reactions:
            model.database.set_reaction(rxnid_ex, reaction_ex)
            reactions.add(rxnid_ex)
        else:
            rxnid_ex = all_reactions[reaction_ex]

        if not model.has_reaction(rxnid_ex):
            added.add(rxnid_ex)
        model.add_reaction(rxnid_ex)
        added_compounds.add(compound)

    return added


def add_all_transport_reactions(model, boundaries, allow_duplicates=False):
    """Add all transport reactions to database and to model.

    Add transport reactions for all boundaries. Boundaries are defined
    by pairs (2-tuples) of compartment IDs. Transport reactions are
    added for all compounds in the model, not just for compounds in the
    two boundary compartments.

    Args:
        model: :class:`psamm.metabolicmodel.MetabolicModel`.
        boundaries: Set of compartment boundary pairs.

    Returns:
        Set of IDs of reactions that were added.
    """

    all_reactions = {}
    if not allow_duplicates:
        # TODO: Avoid adding reactions that already exist in the database.
        # This should be integrated in the database.
        for rxnid in model.database.reactions:
            rx = model.database.get_reaction(rxnid)
            all_reactions[rx] = rxnid

    boundary_pairs = set()
    for source, dest in boundaries:
        if source != dest:
            boundary_pairs.add(tuple(sorted((source, dest))))

    added = set()
    added_pairs = set()
    initial_compounds = set(model.compounds)
    reactions = set(model.database.reactions)
    for compound in initial_compounds:
        for c1, c2 in boundary_pairs:
            compound1 = compound.in_compartment(c1)
            compound2 = compound.in_compartment(c2)
            pair = compound1, compound2
            if pair in added_pairs:
                continue

            rxnid_tp = create_transport_id(reactions, compound1, compound2)

            reaction_tp = Reaction(Direction.Both, {
                compound1: -1,
                compound2: 1
            })
            if reaction_tp not in all_reactions:
                model.database.set_reaction(rxnid_tp, reaction_tp)
                reactions.add(rxnid_tp)
            else:
                rxnid_tp = all_reactions[reaction_tp]

            if not model.has_reaction(rxnid_tp):
                added.add(rxnid_tp)
            model.add_reaction(rxnid_tp)
            added_pairs.add(pair)

    return added


def create_extended_model(model, db_penalty=None, ex_penalty=None,
                          tp_penalty=None, penalties=None):
    """Create an extended model for gap-filling.

    Create a :class:`psamm.metabolicmodel.MetabolicModel` with
    all reactions added (the reaction database in the model is taken
    to be the universal database) and also with artificial exchange
    and transport reactions added. Return the extended
    :class:`psamm.metabolicmodel.MetabolicModel`
    and a weight dictionary for added reactions in that model.

    Args:
        model: :class:`psamm.datasource.native.NativeModel`.
        db_penalty: penalty score for database reactions, default is `None`.
        ex_penalty: penalty score for exchange reactions, default is `None`.
        tb_penalty: penalty score for transport reactions, default is `None`.
        penalties: a dictionary of penalty scores for database reactions.
    """

    # Create metabolic model
    model_extended = model.create_metabolic_model()
    extra_compartment = model.extracellular_compartment

    compartment_ids = set(c.id for c in model.compartments)

    # Add database reactions to extended model
    if len(compartment_ids) > 0:
        logger.info(
            'Using all database reactions in compartments: {}...'.format(
                ', '.join('{}'.format(c) for c in compartment_ids)))
        db_added = add_all_database_reactions(model_extended, compartment_ids)
    else:
        logger.warning(
            'No compartments specified in the model; database reactions will'
            ' not be used! Add compartment specification to model to include'
            ' database reactions for those compartments.')
        db_added = set()

    # Add exchange reactions to extended model
    logger.info(
        'Using artificial exchange reactions for compartment: {}...'.format(
            extra_compartment))
    ex_added = add_all_exchange_reactions(
        model_extended, extra_compartment, allow_duplicates=False)

    # Add transport reactions to extended model
    boundaries = model.compartment_boundaries
    if len(boundaries) > 0:
        logger.info(
            'Using artificial transport reactions for the compartment'
            ' boundaries: {}...'.format(
                '; '.join('{}<->{}'.format(c1, c2) for c1, c2 in boundaries)))
        tp_added = add_all_transport_reactions(
            model_extended, boundaries, allow_duplicates=True)
    else:
        logger.warning(
            'No compartment boundaries specified in the model;'
            ' artificial transport reactions will not be used!')
        tp_added = set()

    # Add penalty weights on reactions
    weights = {}
    if db_penalty is not None:
        weights.update((rxnid, db_penalty) for rxnid in db_added)
    else:
        weights.update((rxnid, 1) for rxnid in db_added)
    if tp_penalty is not None:
        weights.update((rxnid, tp_penalty) for rxnid in tp_added)
    else:
        weights.update((rxnid, 1) for rxnid in tp_added)
    if ex_penalty is not None:
        weights.update((rxnid, ex_penalty) for rxnid in ex_added)
    else:
        weights.update((rxnid, 1) for rxnid in ex_added)

    if penalties is not None:
        for rxnid, penalty in iteritems(penalties):
            weights[rxnid] = penalty
    return model_extended, weights
