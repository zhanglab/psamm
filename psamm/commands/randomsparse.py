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
# Copyright 2015  Keith Dufault-Thompson <keitht547@my.uri.edu>

from __future__ import unicode_literals

import time
import random
import logging

from ..command import Command, MetabolicMixin, SolverCommandMixin, CommandError
from .. import fluxanalysis, util
from ..expression import boolean

from six import string_types

logger = logging.getLogger(__name__)


class RandomSparseNetworkCommand(MetabolicMixin, SolverCommandMixin, Command):
    """Find a random minimal network of model reactions.

    Given a reaction to optimize and a threshold, delete reactions randomly
    until the flux of the reaction to optimize falls under the threshold.
    Keep deleting reactions until no more reactions can be deleted. By default
    this uses standard FBA (not tFBA). Since the internal fluxes are irrelevant
    the FBA and tFBA are equivalent for this purpose.

    The threshold can be specified as an absolute flux (e.g. '1.23') or a
    relative flux of the full model flux (e.g. '40.5%').
    """

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--tfba', help='Enable thermodynamic constraints on FBA',
            action='store_true', default=False)
        parser.add_argument('--objective', help='Reaction flux to maximize')
        parser.add_argument(
            'threshold', help='Threshold of max reaction flux',
            type=util.MaybeRelative)
        parser.add_argument(
            '--type', help='Type of deletion to perform',
            choices=['reactions', 'exchange', 'genes'], type=str,
            required=True)
        super(RandomSparseNetworkCommand, cls).init_parser(parser)

    def run(self):
        if self._args.objective is not None:
            reaction = self._args.objective
        else:
            reaction = self._model.get_biomass_reaction()
            if reaction is None:
                raise CommandError('The biomass reaction was not specified')

        if not self._mm.has_reaction(reaction):
            raise CommandError('Specified reaction is not in model: {}'.format(
                reaction))

        if not self._args.tfba:
            solver = self._get_solver()
        else:
            solver = self._get_solver(integer=True)

        p = fluxanalysis.FluxBalanceProblem(self._mm, solver)

        if self._args.tfba:
            p.add_thermodynamic()

        threshold = self._args.threshold
        if threshold.relative:
            p.maximize(reaction)
            threshold.reference = p.get_flux(reaction)

        flux_threshold = float(threshold)

        logger.info('Flux threshold for {} is {}'.format(
            reaction, flux_threshold))

        if self._args.type in ('reactions', 'exchange'):
            self._delete_reactions(p, reaction, flux_threshold)
        elif self._args.type == 'genes':
            self._delete_genes(p, reaction, flux_threshold)

    def _delete_reactions(self, prob, obj_reaction, flux_threshold):
        essential = {obj_reaction}
        if self._args.type == 'exchange':
            reactions = set()
            for reaction_id in self._mm.reactions:
                if self._mm.is_exchange(reaction_id):
                    reactions.add(reaction_id)
        else:
            reactions = set(self._mm.reactions)
        test_set = reactions - essential
        deleted = set()

        start_time = time.time()

        while len(test_set) > 0:
            testing_reaction = random.sample(test_set, 1)[0]
            test_set.remove(testing_reaction)

            flux_var = prob.get_flux_var(testing_reaction)
            c, = prob.prob.add_linear_constraints(flux_var == 0)

            logger.info('Trying FBA without reaction {}...'.format(
                testing_reaction))

            try:
                prob.maximize(obj_reaction)
            except fluxanalysis.FluxBalanceError:
                logger.info(
                    'FBA is infeasible, marking {} as essential'.format(
                        testing_reaction))
                c.delete()
                essential.add(testing_reaction)
                continue

            logger.debug('Reaction {} has flux {}'.format(
                obj_reaction, prob.get_flux(obj_reaction)))

            if prob.get_flux(obj_reaction) < flux_threshold:
                c.delete()
                essential.add(testing_reaction)
                logger.info('Reaction {} was essential'.format(
                    testing_reaction))
            else:
                deleted.add(testing_reaction)
                logger.info('Reaction {} was deleted'.format(testing_reaction))

        logger.info('Solving took {:.2f} seconds'.format(
            time.time() - start_time))

        for reaction_id in sorted(reactions):
            value = 0 if reaction_id in deleted else 1
            print('{}\t{}'.format(reaction_id, value))

        logger.info('Essential reactions: {}/{}'.format(
            len(essential), len(reactions)))

    def _delete_genes(self, prob, obj_reaction, flux_threshold):
        genes = set()
        gene_assoc = {}
        for reaction in self._model.parse_reactions():
            assoc = None
            if reaction.genes is None:
                continue
            elif isinstance(reaction.genes, string_types):
                assoc = boolean.Expression(reaction.genes)
            else:
                variables = [boolean.Variable(g) for g in reaction.genes]
                assoc = boolean.Expression(boolean.And(*variables))
            genes.update(v.symbol for v in assoc.variables)
            gene_assoc[reaction.id] = assoc

        essential = set()
        deleted = set()
        test_set = set(genes)
        reactions = set(self._mm.reactions)

        start_time = time.time()

        while len(test_set) > 0:
            testing_gene = random.sample(test_set, 1)[0]
            test_set.remove(testing_gene)
            new_gene_assoc = {}
            deleted_reactions = set()

            logger.info('Trying model without gene {}...'.format(
                testing_gene))

            for reaction in reactions:
                if reaction not in gene_assoc:
                    continue
                assoc = gene_assoc[reaction]
                if boolean.Variable(testing_gene) in assoc.variables:
                    new_assoc = assoc.substitute(
                        lambda v: v if v.symbol != testing_gene else False)
                    if new_assoc.has_value() and not new_assoc.value:
                        deleted_reactions.add(reaction)
                    else:
                        new_gene_assoc[reaction] = new_assoc
                else:
                    new_gene_assoc[reaction] = assoc

            if obj_reaction in deleted_reactions:
                logger.info(
                    'Marking gene {} as essential because the objective'
                    ' reaction depends on this gene...'.format(testing_gene))
                essential.add(testing_gene)
                continue

            if len(deleted_reactions) == 0:
                logger.info(
                    'No reactions were removed when gene {}'
                    ' was deleted'.format(testing_gene))
                deleted.add(testing_gene)
                gene_assoc = new_gene_assoc
                continue

            logger.info(
                'Reactions removed when gene {} is deleted: {}'.format(
                    testing_gene, len(deleted_reactions)))
            logger.info('Deleted reactions: {}'.format(
                ', '.join(deleted_reactions)))

            constraints = []
            for reaction in deleted_reactions:
                flux_var = prob.get_flux_var(reaction)
                c, = prob.prob.add_linear_constraints(flux_var == 0)
                constraints.append(c)

            try:
                prob.maximize(obj_reaction)
            except fluxanalysis.FluxBalanceError:
                logger.info(
                    'FBA is infeasible, marking {} as essential'.format(
                        testing_gene))
                for c in constraints:
                    c.delete()
                essential.add(testing_gene)
                continue

            logger.debug('Reaction {} has flux {}'.format(
                obj_reaction, prob.get_flux(obj_reaction)))

            if prob.get_flux(obj_reaction) < flux_threshold:
                for c in constraints:
                    c.delete()
                essential.add(testing_gene)
                logger.info('Gene {} was essential'.format(testing_gene))
            else:
                deleted.add(testing_gene)
                reactions.difference_update(deleted_reactions)
                gene_assoc = new_gene_assoc
                logger.info('Gene {} was deleted'.format(testing_gene))

        logger.info('Solving took {:.2f} seconds'.format(
            time.time() - start_time))

        for gene_id in sorted(genes):
            value = 0 if gene_id in deleted else 1
            print('{}\t{}'.format(gene_id, value))

        logger.info('Essential genes: {}/{}'.format(
            len(essential), len(genes)))

        logger.info('Reactions in minimal network: {}'.format(
            ', '.join(sorted(reactions))))
