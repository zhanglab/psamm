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
import logging

from ..command import (Command, MetabolicMixin, ObjectiveMixin,
                       SolverCommandMixin, FilePrefixAppendAction)
from .. import fluxanalysis
from .. import moma
from ..expression import boolean

from six import string_types

logger = logging.getLogger(__name__)


class GeneDeletionCommand(MetabolicMixin, ObjectiveMixin, SolverCommandMixin,
                          Command):
    """Find reactions requiring a specified gene and delete them from model.

    Reports the new objective flux after the gene was deleted.
    """
    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--gene', metavar='genes', action=FilePrefixAppendAction,
            type=str, default=[], help='Delete multiple genes from model')
        parser.add_argument(
            '--method', metavar='method', choices=['fba', 'lp2', 'lp3', 'qlp2', 'qlp3'],
            type=str, default='fba', help='Select which method to solve the problem with')
        super(GeneDeletionCommand, cls).init_parser(parser)

    def run(self):
        obj_reaction = self._get_objective()

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

        reactions = set(self._mm.reactions)
        start_time = time.time()
        testing_genes = set(self._args.gene)
        deleted_reactions = set()

        logger.info('Trying model without genes: {}...'.format(
                    ', '.join(sorted(testing_genes))))

        for reaction in reactions:
            if reaction not in gene_assoc:
                continue
            assoc = gene_assoc[reaction]
            if any(boolean.Variable(gene) in assoc.variables
                    for gene in testing_genes):
                new_assoc = assoc.substitute(
                    lambda v: v if v.symbol not in testing_genes else False)
                if new_assoc.has_value() and not new_assoc.value:
                    logger.info('Deletion reaction {}...'.format(reaction))
                    deleted_reactions.add(reaction)

        solver = self._get_solver()

        # Solve with FBA
        if self._args.method == 'fba':
            logger.info('Solving using FBA...')
            # Set up the FBA problem
            prob = fluxanalysis.FluxBalanceProblem(self._mm, solver)

            # Try to maximize the biomass
            try:
                prob.maximize(obj_reaction)
            except fluxanalysis.FluxBalanceError as e:
                self.report_flux_balance_error(e)

            # Get the wildtype flux
            wild = prob.get_flux(obj_reaction)

            # Add linear constraints to the reactions for the genes that the user specifies
            for reaction in deleted_reactions:
                flux_var = prob.get_flux_var(reaction)
                prob.prob.add_linear_constraints(flux_var == 0)

            # Maximize the biomass for the wild type
            prob.maximize(obj_reaction)
            deleteflux = prob.get_flux(obj_reaction)


            if wild != 0:
                logger.info(
                    'Objective reaction has {:.2%} flux of wild type flux'.format(
                        abs(deleteflux / wild)))
                logger.info(
                    'Solving took {:.2f} seconds'.format(time.time() - start_time))
                logger.info(
                    'Objective reaction after gene deletion has flux {}'.format(
                        abs(deleteflux) if deleteflux == 0 else deleteflux))

        # Solve with MOMA LP3
        elif self._args.method in ['lp2', 'lp3', 'qlp2', 'qlp3']:

            # Set up the MOMA solver
            p = moma.MOMAProblem(self._mm, solver)

            # Get the wildtype biomass for comparison later
            wild = p.get_fba_biomass(obj_reaction)

            # Get the wildtype fluxes that we are going to constrain are model with
            wt_fluxes = p.get_fba_flux(obj_reaction)
            #print(deleted_reactions)
            constr = []

            # Set all the reactions assiciated to the gene to a flux state of 0
            for reaction in deleted_reactions:
                constr.extend(p.prob.add_linear_constraints(
                    p.get_flux_var(reaction) == 0))

            if self._args.method == 'lp3':
                logger.info('Solving using MOMA LP3...')
                p.minimize_l1(obj_reaction, wild, wt_fluxes)
            elif self._args.method == 'lp2':
                logger.info('Solving using MOMA LP2...')
                p.minimize_lp2(obj_reaction, wt_fluxes)
            elif self._args.method == 'qlp2':
                logger.info('Solving using MOMA QLP2...')
                p.minimize_qlp2(obj_reaction, wt_fluxes)
            elif self._args.method == 'qlp3':
                logger.info('Solving using MOMA QLP3...')
                p.test_minimize_l2(obj_reaction)

            try:
                deleteflux = p.get_flux(obj_reaction)
                logger.info(
                    'Objective reaction has {:.2%} flux of wild type flux'.format(
                        abs(deleteflux / wild)))
                logger.info(
                    'Solving took {:.2f} seconds'.format(time.time() - start_time))
                logger.info(
                    'Objective reaction after gene deletion has flux {}'.format(
                        abs(deleteflux) if deleteflux == 0 else deleteflux))
            except moma.MOMAError:
                logger.info('Error computing the flux')

            # Do upon exit
            finally:
                for c in constr:
                    c.delete()
