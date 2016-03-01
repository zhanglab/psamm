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
# Copyright 2016  Keith Dufault-Thompson <keitht547@my.uri.edu>


from __future__ import unicode_literals
import time
import logging

from ..command import SolverCommandMixin, MetabolicMixin, Command, CommandError
from .. import fluxanalysis
from psamm.datasource.native import NativeModel
import csv

# = logging.getLogger(__name__)


class MultimediaCommand(MetabolicMixin, SolverCommandMixin, Command):
    """Export the metabolic model as an Excel workbook"""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--carbon', type=str, help='File path to carbon source id list')
        parser.add_argument(
            '--electron', type=str, help='File path to electron acceptor source id list')
        parser.add_argument(
            '--loop-removal', help='Select type of loop removal constraints',
            choices=['none', 'tfba', 'l1min'], default='none')
        parser.add_argument(
            '--all-reactions', help='Show all reaction fluxes',
            action='store_true')
        parser.add_argument(
            '--epsilon', type=float, help='Threshold for flux minimization',
            default=1e-5)
        parser.add_argument('reaction', help='Reaction to maximize', nargs='?')
        super(MultimediaCommand, cls).init_parser(parser)

    def run(self):
        carbon_list = open(self._args.carbon, mode='rU')
        electron_list = open(self._args.electron, mode='rU')

        carbon_l = []
        electron_l = []
        for carbon in carbon_list.readlines():
            carbon = carbon.rstrip()
            carbon_l.append(carbon)
        for electron in electron_list.readlines():
            electron = electron.rstrip()
            electron_l.append(electron)



        compound_name = {}
        for compound in self._model.parse_compounds():
            if 'name' in compound.properties:
                compound_name[compound.id] = compound.properties['name']
            elif compound.id not in compound_name:
                compound_name[compound.id] = compound.id

        # Reaction genes information
        reaction_genes = {}
        for reaction in self._model.parse_reactions():
            if 'genes' in reaction.properties:
                reaction_genes[reaction.id] = reaction.properties['genes']

        if self._args.reaction is not None:
            reaction = self._args.reaction
        else:
            reaction = self._model.get_biomass_reaction()
            if reaction is None:
                raise CommandError('The biomass reaction was not specified')

        if not self._mm.has_reaction(reaction):
            raise CommandError('Specified reaction is not in model: {}'.format(
                reaction))

        #logger.info('Using {} as objective'.format(reaction))

        for csource in carbon_l:
            for eaccept in electron_l:
                output_file = open('{}.tsv'.format(eaccept), mode='a')

                result = self.run_carbon_fba(reaction, csource, eaccept)
                #output_file.write('{}\t{}\t{}\t{}\n'.format('', 'Carbon_flux', 'Electron_acc', 'Core_Biomass'))
                for c_e, c_l, c_u, e_e, e_l, e_u, b_e, b_l, b_u in result:
                    output_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(csource, c_e, c_l, c_u,
                                                                                  e_e, e_l, e_u, b_e, b_l, b_u))





                '''
                optimum = None
                total_reactions = 0
                nonzero_reactions = 0
                output_file = open('{}_{}.tsv'.format(eaccept, csource), mode='w')
                for reaction_id, flux in sorted(result):
                    total_reactions += 1
                    if abs(flux) > self._args.epsilon:
                        nonzero_reactions += 1

                    if abs(flux) > self._args.epsilon or self._args.all_reactions:
                        rx = self._mm.get_reaction(reaction_id)
                        rx_trans = rx.translated_compounds(
                            lambda x: compound_name.get(x, x))
                        genes = reaction_genes.get(reaction_id, '')

                        output_file.write('{}\t{}\t{}\t{}\n'.format(
                            reaction_id, flux, rx_trans, genes))


                    # Remember flux of requested reaction
                    if reaction_id == reaction:
                        optimum = flux

                logger.info('Objective flux: {}'.format(optimum))
                logger.info('Reactions at zero flux: {}/{}'.format(
                    total_reactions - nonzero_reactions, total_reactions))
                '''


    def run_carbon_fba(self, reaction, carbon, electron):
        """Run standard FBA on model."""
        solver = self._get_solver()
        model = self._mm
        model_c = model.copy()
        model_c.limits['EX_{}_e'.format(carbon)]._assign_lower(-10)
        model_c.limits['EX_{}_e'.format(electron)]._assign_lower(-1000)
        if electron == 'cpd_fe3':
            model_c.limits['EX_cpd_fe2_e']._assign_upper(1000)
        p = fluxanalysis.FluxBalanceProblem(model_c, solver)
        #start_time = time.time()

        p.maximize(reaction)
        obj_flux = p.get_flux(reaction)
        c_list = []
        e_list = []
        b_list = []
        c_list.append('EX_{}_e'.format(carbon))
        e_list.append('EX_{}_e'.format(electron))
        b_list.append(reaction)
        for exchange_c, (lower_c, upper_c) in fluxanalysis.flux_variability(model_c, c_list, {reaction: obj_flux}, tfba=False, solver=self._get_solver()):
            c_e = exchange_c
            c_l = lower_c
            c_u = upper_c
        for exchange_e, (lower_e, upper_e) in fluxanalysis.flux_variability(model_c, e_list, {reaction: obj_flux}, tfba=False, solver=self._get_solver()):
            e_e = exchange_e
            e_l = lower_e
            e_u = upper_e
        for exchange_b, (lower_b, upper_b) in fluxanalysis.flux_variability(model_c, b_list, {reaction: obj_flux}, tfba=False, solver=self._get_solver()):
            b_e = exchange_b
            b_l = lower_b
            b_u = upper_b
        yield(c_e, c_l, c_u, e_e, e_l, e_u, b_e, b_l, b_u)
        #logger.info('Solving took {:.2f} seconds'.format(
        #    time.time() - start_time))

        #yield 'EX_{}_e'.format(carbon), p.get_flux_var('EX_{}_e'.format(carbon)), 'EX_{}_e'.format(electron), p.get_flux_var('EX_{}_e'.format(electron)), reaction, p.get_flux_var(reaction)
        #for reaction_id in self._mm.reactions:
        #    yield reaction_id, p.get_flux(reaction_id)

