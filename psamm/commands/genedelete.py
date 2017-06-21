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
            '--method', metavar='method',
            choices=['fba', 'lin_moma', 'lin_moma2', 'moma', 'moma2'],
            type=str, default='fba',
            help='Select which method to use. (fba, lin_moma, lin_moma2, '
                 'moma, moma2)')
        super(GeneDeletionCommand, cls).init_parser(parser)

    def run(self):
        """Delete the specified gene and solve using the desired method."""
        obj_reaction = self._get_objective()

        genes = set()
        gene_assoc = {}
        for reaction in self._model.reactions:
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
                    logger.info('Deleting reaction {}...'.format(reaction))
                    deleted_reactions.add(reaction)

        if self._args.method in ['moma', 'moma2']:
            solver = self._get_solver(quadratic=True)
        else:
            solver = self._get_solver()

        if self._args.method == 'fba':
            logger.info('Solving using FBA...')
            prob = fluxanalysis.FluxBalanceProblem(self._mm, solver)

            try:
                prob.maximize(obj_reaction)
            except fluxanalysis.FluxBalanceError as e:
                self.report_flux_balance_error(e)

            wild = prob.get_flux(obj_reaction)

            for reaction in deleted_reactions:
                flux_var = prob.get_flux_var(reaction)
                prob.prob.add_linear_constraints(flux_var == 0)

            prob.maximize(obj_reaction)
            deleteflux = prob.get_flux(obj_reaction)

            if wild != 0:
                logger.info(
                    'Objective reaction has {:.2%} '
                    'flux of wild type flux'.format(abs(deleteflux / wild)))
                logger.info(
                    'Solving took {:.2f} seconds'.format(
                        time.time() - start_time))
                logger.info(
                    'Objective reaction after gene deletion '
                    'has flux {}'.format(
                        abs(deleteflux) if deleteflux == 0 else deleteflux))

        elif self._args.method in ['lin_moma', 'lin_moma2', 'moma', 'moma2']:
            p = moma.MOMAProblem(self._mm, solver)
            wild = p.get_fba_obj_flux(obj_reaction)
            wt_fluxes = p.get_minimal_fba_flux(obj_reaction)

            constr = []
            try:
                for reaction in deleted_reactions:
                    constr.extend(p.prob.add_linear_constraints(
                        p.get_flux_var(reaction) == 0))

                if self._args.method == 'lin_moma2':
                    logger.info('Solving using linear moma2...')
                    p.lin_moma2(obj_reaction, wild)
                elif self._args.method == 'lin_moma':
                    logger.info('Solving using linear moma...')
                    p.lin_moma(wt_fluxes)
                elif self._args.method == 'moma':
                    logger.info('Solving using moma...')
                    p.moma(wt_fluxes)
                elif self._args.method == 'moma2':
                    logger.info('Solving using moma2...')
                    p.moma2(obj_reaction, wild)

                deleteflux = p.get_flux(obj_reaction)
                logger.info(
                    '''Objective reaction has {:.2%} flux of wild
                    type flux'''.format(abs(deleteflux / wild)))
                logger.info(
                    'Solving took {:.2f} seconds'.format(
                        time.time() - start_time))
                logger.info(
                    'Objective reaction after gene deletion '
                    'has flux {}'.format(
                        abs(deleteflux) if deleteflux == 0 else deleteflux))

            except moma.MOMAError:
                logger.info('Error computing the flux')
            # Reset the model for the next gene
            finally:
                for c in constr:
                    c.delete()
