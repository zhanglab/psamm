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
# Copyright 2016  Keith Dufault-Thompson <keitht547@my.uri.edu>
# Copyright 2016  Julie Cuddigan <julie_cuddigan@my.uri.edu>

"""Metabolic Adjustment by Differential Expression (MADE) command."""

from __future__ import unicode_literals

import time
import logging
import math
import csv

from six import iteritems

from ..command import SolverCommandMixin, MetabolicMixin, Command
from ..util import MaybeRelative
from ..lpsolver import lp
from ..expression import boolean

# Module-level logging
logger = logging.getLogger(__name__)


class MadeFluxBalance(MetabolicMixin, SolverCommandMixin, Command):
    """Run MADE flux balance analysis on the model.

    Args:
        gene_var1 = Dictionary, key:value = gene expression objects:their new
            variable id, first set
        gene_var2 = Dictionary; key:value = gene expression objects:their new
            variable id, second set
        var_ineqvar1 = xi; Dictionary, key:value = new variable ids:their
            defined inequality variable, first set
        var_ineqvar2 = xi+1; Dictionary, key:value = new variable ids:their
            defined inequality variable, second set
        gene_pval = Dictionary, key:value = gene ID:gene fold change
            probability (pvalue)
        gene_diff = Dictionary, key:value = gene ID: binary up/down/constant
            regulation values
        gvdict = Dictionary, key:value = gene ID:defined variable ids from both
            sets (each key has 2 values)
        problem = Flux balance problem
    """

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--flux-threshold',
            help='Enter maximum objective flux as a decimal or percent',
            type=MaybeRelative, default=MaybeRelative('100%'))
        parser.add_argument(
            '--transc-file', help='Enter path to transcriptomic data file',
            metavar='FILE')
        parser.add_argument('--fva', help='Enable FVA', action='store_true')
        super(MadeFluxBalance, cls).init_parser(parser)

    def run(self):
        """Run MADE implementation."""
        gene_dict = self.get_gene_dict()

        biomass_fun = self._model.biomass_reaction

        # Create problem instance
        solver = self._get_solver(integer=True)
        prob = solver.create_problem()
        v_0 = prob.namespace()
        v_1 = prob.namespace()

        # Define flux variables
        for reaction_id in self._mm.reactions:
            lower, upper = self._mm.limits[reaction_id]
            v_0.define([reaction_id], lower=lower, upper=upper)
            v_1.define([reaction_id], lower=lower, upper=upper)

        # Create mass balance constraints for both conditions
        massbalance_0_lhs = {compound: 0 for compound in self._mm.compounds}
        massbalance_1_lhs = {compound: 0 for compound in self._mm.compounds}
        for (compound, reaction_id), value in iteritems(self._mm.matrix):
            massbalance_0_lhs[compound] += v_0(reaction_id) * value
            massbalance_1_lhs[compound] += v_1(reaction_id) * value
        for _, lhs in iteritems(massbalance_0_lhs):
            prob.add_linear_constraints(lhs == 0)
        for _, lhs in iteritems(massbalance_1_lhs):
            prob.add_linear_constraints(lhs == 0)

        start_time = time.time()

        # Set biomass flux threshold
        flux_threshold = self._args.flux_threshold
        if flux_threshold.relative:
            prob.set_objective(v_0(biomass_fun))
            result = prob.solve(lp.ObjectiveSense.Maximize)
            if not result:
                raise Exception('Failed to solve FBA')
            flux_threshold.reference = result.get_value(v_0(biomass_fun))

        prob.add_linear_constraints(v_0(biomass_fun) >= float(flux_threshold))
        prob.add_linear_constraints(v_1(biomass_fun) >= float(flux_threshold))

        gene_term_0 = prob.namespace()
        gene_term_1 = prob.namespace()

        reaction_0 = prob.namespace(
            self._mm.reactions, types=lp.VariableType.Binary)
        reaction_1 = prob.namespace(
            self._mm.reactions, types=lp.VariableType.Binary)

        for rxn_id, exp in sorted(iteritems(gene_dict)):
            create_gpr_constraints(
                prob, rxn_id, exp, reaction_0(rxn_id), gene_term_0)
            create_gpr_constraints(
                prob, rxn_id, exp, reaction_1(rxn_id), gene_term_1)

        if self._args.transc_file is not None:
            con1, con2, gene_pval, gene_diff = idc(
                open_file(self._args.transc_file))

        add_final_constraints(self._mm, prob, v_0, reaction_0)
        add_final_constraints(self._mm, prob, v_1, reaction_1)
        result = make_obj_fun(
            prob, gene_diff, gene_pval, gene_term_0, gene_term_1)

        # Run FBA
        for reaction_id in sorted(self._mm.reactions):
            print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                reaction_id, result.get_value(v_0(reaction_id)),
                result.get_value(v_1(reaction_id)),
                result.get_value(reaction_0(reaction_id)) > 0.5,
                result.get_value(reaction_1(reaction_id)) > 0.5,
                self._mm.get_reaction(reaction_id),
                gene_dict.get(reaction_id, '')))

        logger.info('Solving took {:.2f} seconds'.format(
            time.time() - start_time))

    def get_gene_dict(self):
        """Using the reaction file called inside of the model file, it returns
        a dictionary with reaction IDs as keys and their associated
        gene-protein reaction (GPR) logic (i.e. (gene 1 and gene 2) or gene 3)
        as values of type str.
        """
        gene_dict = {}
        for reaction in self._model.parse_reactions():
            if reaction.genes is not None:
                gene_dict[reaction.id] = boolean.Expression(reaction.genes)

        return gene_dict


def make_obj_fun(prob, gene_diff, gene_pval, gene_term_0, gene_term_1):
    """Constructs the MADE objective funtion from dictionaries of LP variables.

    Objective function consists of the summation of three functions dependent
    on the up/down regulation of gene expression between conditions. The
    functions contain a weighting function, and the difference between the
    binary representations of condition 1 and condition 2.
    """
    i_vars = 0.0  # Increasing gene expression
    d_vars = 0.0  # Decreasing gene expression
    c_vars = 0.0  # Constant gene expression

    def weight(p):
        return -math.log10(p)

    for gene in gene_pval:
        # Comment by Julie/Matt?: Limitation of math.log()
        wp = max(2.2204460492e-16, gene_pval[gene])

        x_0 = gene_term_0(boolean.Variable(gene))
        x_1 = gene_term_1(boolean.Variable(gene))
        # x_delta = x_0 XOR X_1
        prob.define(('xor', gene), types=lp.VariableType.Binary)
        x_delta = prob.var(('xor', gene))
        prob.add_linear_constraints(
            x_delta <= x_0 + x_1,
            x_delta >= x_0 - x_1, x_delta >= x_1 - x_0,
            x_delta <= 2 - x_0 - x_1)

        if gene_diff[gene] == 1:
            i_vars += weight(wp) * (x_1 - x_0)
        elif gene_diff[gene] == -1:
            d_vars += weight(wp) * (x_0 - x_1)
        elif gene_diff[gene] == 0:
            c_vars += weight(wp) * x_delta

    objective = i_vars + d_vars - c_vars

    prob.set_objective(objective)
    result = prob.solve(lp.ObjectiveSense.Maximize)
    if not result:
        raise Exception('Unable to solve: ' + result.status)

    obj_value = result.get_value(objective)
    logger.info('Objective: {}'.format(obj_value))

    return result


def create_gpr_constraints(prob, rxn_id, exp, reaction_var, gene_term):
    """Opens all gene-logic containers, defines content, outputs the linear
    inequalities by calling bool_ineqs().  Sorts data into dictionaries that
    are used in other functions.  Is recursive. No output.

    Args:
        exp_obj: All of the expression objects (genes, AND, OR)
        var_gen: Counter used for relabeling the genes and arguments as
            variables
        new_var_id: Variable ID, also includes original reaction ID for first
            layer
    """
    next_terms = iter([exp.root])
    current_type = None
    variable = reaction_var
    arguments = []
    stack = []
    while True:
        try:
            term = next(next_terms)
        except StopIteration:
            term = None

        if gene_term.has_variable(term):
            arguments.append(gene_term(term))
        elif isinstance(term, boolean.Variable):
            gene_term.define([term], types=lp.VariableType.Binary)
            term_var = gene_term(term)
            arguments.append(term_var)
        elif isinstance(term, (boolean.And, boolean.Or)):
            stack.append((current_type, next_terms, variable, arguments))
            current_type = term.__class__
            next_terms = iter(term)
            arguments = []
            gene_term.define([term], types=lp.VariableType.Binary)
            variable = gene_term(term)
        else:
            # End of term
            if current_type is None:
                prob.add_linear_constraints(variable == arguments[0])
                break
            elif current_type == boolean.And:
                prob.add_linear_constraints(
                    *and_constraints(variable, arguments))
            elif current_type == boolean.Or:
                prob.add_linear_constraints(
                    *or_constraints(variable, arguments))

            term_var = variable
            current_type, next_terms, variable, arguments = stack.pop()
            arguments.append(term_var)


def and_constraints(var, arguments):
    """Create constraints for boolean AND.

    Creates constraints for: var = And(arguments) where var and each argument
    is a binary variable.
    """
    var_sum = 0
    for arg in arguments:
        yield var <= arg
        var_sum += arg
    yield var >= var_sum - (len(arguments) - 1)


def or_constraints(var, arguments):
    """Create constraints for boolean OR.

    Creates constraints for: var = Or(arguments) where var and each argument
    is a binary variable.
    """
    var_sum = 0
    for arg in arguments:
        yield var >= arg
        var_sum += arg
    yield var <= var_sum


def open_file(path):
    """Returns the contents of model file in a tuple of dictionaries.
    File Form: tsv format, FOUR Columns: (1) Gene name, (2) Condition 1 Data,
    (3) Condition 2 Data, (4) P-value of the fold change for transition 1->2.
    """
    file1 = open(path)
    con1_dict = {}
    con2_dict = {}
    pval_dict = {}

    file1.readline()
    for row in csv.reader(file1, delimiter=str('\t')):
        con1_dict[row[0]] = float(row[1])
        con2_dict[row[0]] = float(row[2])
        if float(row[3]) == float(0.0):
            pval_dict[row[0]] = 1e-400
        else:
            pval_dict[row[0]] = float(row[3])

    return con1_dict, con2_dict, pval_dict


def idc(dicts):
    """Used for accessing the list of dictionaries created in open_file()
    Creates a dictionary for the gene ID and a value = [-1, 0, +1]
    corresponding to decreasing, constant, and inreasing expression between the
    conditions.
    """
    con1 = dicts[0]
    con2 = dicts[1]
    pval = dicts[2]
    diff = {}

    for key in con1:
        if con2[key] == con1[key]:
            diff[key] = 0
        elif con2[key] > con1[key]:
            diff[key] = 1
        else:
            diff[key] = -1

    return con1, con2, pval, diff


def add_final_constraints(mm, prob, v, z):
    """Takes the metabolic model, the LP Problem, and the binary
    dictionaries of each condition.  Adds constraints connecting flux
    variables, reactions, and their flux bounds.
    """
    for rxn in mm.reactions:
        vmin = mm.limits[rxn].lower
        vmax = mm.limits[rxn].upper
        flux_var = v(rxn)
        active = z(rxn)

        prob.add_linear_constraints(
            active*vmax >= flux_var, flux_var >= active*vmin)
