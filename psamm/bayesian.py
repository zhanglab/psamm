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

"""Calculate model mapping likelihood with bayesian."""
from __future__ import print_function, division
from future.utils import itervalues

from builtins import object, range
from itertools import product
from multiprocessing import Pool
import operator
import sys
from collections import namedtuple, defaultdict
import pandas as pd
import time

import psamm.bayesian_util as util
from psamm.formula import Formula
from psamm.expression.boolean import Expression
from functools import reduce

import logging
logger = logging.getLogger(__name__)


CompoundEntry = namedtuple(
    'CompoundEntry',
    ['id', 'name', 'formula', 'charge', 'kegg', 'cas'])

ReactionEntry = namedtuple(
    'ReactionEntry',
    ['id', 'name', 'genes', 'equation', 'subsystem', 'ec'])


class MappingModel(object):
    """Generate the internal structure for model mapping.

    Args:
        model: :class:`psamm.datasource.native.NativeModel`
    """

    def __init__(self, model):
        self._name = model.name
        self._read_compounds(model.compounds)
        self._read_reactions(model.reactions)

        self._genes = set()
        self._has_gene = list()
        for r in itervalues(self._reactions):
            if r.genes is not None:
                self._has_gene.append(r.id)
                e = Expression(r.genes)
                self._genes.update([g.symbol for g in e.variables])

    def _read_compounds(self, it):
        self._compounds = {}
        for compound in it:
            formula = getattr(compound, 'formula', None)
            if formula is not None:
                formula = Formula.parse(formula)
            else:
                formula = None

            if 'kegg' in compound.properties:
                compound_kegg = compound.properties['kegg']
            else:
                compound_kegg = None

            compound_id = compound.id
            entry = CompoundEntry(
                id=compound_id, name=getattr(compound, 'name', None),
                formula=formula, charge=getattr(compound, 'charge', None),
                kegg=compound_kegg,
                cas=getattr(compound, 'cas', None))
            self._compounds[compound_id] = entry

    def _read_reactions(self, it):
        self._reactions = {}
        for reaction in it:
            genes = getattr(reaction, 'genes', None)
            reaction_id = reaction.id
            equation = getattr(reaction, 'equation', None)
            entry = ReactionEntry(
                id=reaction_id, name=getattr(reaction, 'name', None),
                genes=genes, equation=equation,
                subsystem=getattr(reaction, 'subsystem', None),
                ec=getattr(reaction, 'ec', None))
            self._reactions[reaction_id] = entry

    @property
    def name(self):
        return self._name

    @property
    def reactions(self):
        return self._reactions

    @property
    def compounds(self):
        return self._compounds

    @property
    def genes(self):
        return self._genes

    def print_summary(self):
        """Print model summary"""
        print('Model: {}'.format(self.name))
        print('- Compounds: {}'.format(len(self.compounds)))
        print('- Reactions: {}'.format(len(self.reactions)))
        print('- Genes: {}'.format(len(self.genes)))
        print('- Reactions with gene association: {}'.format(
            len(self._has_gene)))

    def check_reaction_compounds(self):
        """Check that reaction compounds are defined in the model"""
        print('Checking model: {}'.format(self.name))
        consistent = True
        undefined_cpds = []
        for reaction in itervalues(self.reactions):
            if reaction.equation is not None:
                for compound, value in reaction.equation.compounds:
                    if compound.name not in self.compounds:
                        undefined_cpds.append(compound.name)
                        logger.error((
                            '{} in reaction {} is not defined in compound '
                            'list (such as compounds.yaml), please define it '
                            'in compound list before running modelmapping'
                        ).format(compound.name, reaction.id))
                        consistent = False
        return consistent


class BayesianCompoundPredictor(object):
    """Predict model compound mappings based on a Bayesian model.

    Args:
        model1: :class:`psamm.bayesian.MappingModel`.
        model2: :class:`psamm.bayesian.MappingModel`.
        nproc: number of processes used for mapping.
        outpath: the path to output the detailed log file.
        log: whether to output log file of the p_match and p_no_match.
        kegg: whether to compare the KEGG id.
    """

    def __init__(self, model1, model2, nproc=1,
                 outpath='.', log=False, kegg=False):
        self._model1 = model1
        self._model2 = model2
        self._column_list = [
            'p', 'p_id', 'p_name', 'p_charge', 'p_formula', 'p_kegg']
        self._compound_map_p = map_model_compounds(
            self._model1, self._model2, nproc, outpath, log=log, kegg=kegg)

    @property
    def model1(self):
        return self._model1

    @property
    def model2(self):
        return self._model2

    def map(self, c1, c2):
        return self._compound_map_p[0][c1, c2]

    def get_raw_map(self):
        """Return pandas.DataFrame style of raw mapping table."""
        compound_result = pd.DataFrame({
            self._column_list[i]: self._compound_map_p[i]
            for i in range(len(self._column_list))
        })
        return compound_result

    def get_best_map(self, threshold_compound=0):
        """Return :class:`pandas.DataFrame` style of best mapping for each
        query."""
        raw_map = self.get_raw_map()
        query = [i for i, j in raw_map.index.values]
        # use rank instead of idxmax to output multiple top-hitts
        best_index = raw_map.iloc[:, 0].groupby(query).rank(
            method='min', ascending=False).loc[lambda x: x == 1].index.values
        compound_best = raw_map.loc[best_index]
        compound_best.query(
            'p >= @threshold_compound', inplace=True)
        return compound_best

    def get_cpd_pred(self, threshold_compound=0):
        """Return the cpd_pred used for reaction mapping."""
        return self.get_best_map(threshold_compound).loc[:, 'p']


class BayesianReactionPredictor(object):
    """Predict model reaction mappings based on a Bayesian model.

    Args:
        model1: :class:`psamm.bayesian.MappingModel`.
        model2: :class:`psamm.bayesian.MappingModel`.
        cpd_pred: :class:`pandas.Series` with compound pairs as index,
                  compound mapping score as value.
        nproc: number of processes used for mapping.
        outpath: the path to output the detailed log file.
        log: whether to output log file of the p_match and p_no_match.
        gene: whether to compare the gene association.
        compartment_map: dictionary mapping compartment id in the query model
                         to the id in the target model.
        gene_map: dictionary mapping gene id in the query model to the id in
                  the target model.
    """

    def __init__(self, model1, model2, cpd_pred, nproc=1,
                 outpath='.', log=False, gene=False,
                 compartment_map={}, gene_map={}):
        self._model1 = model1
        self._model2 = model2
        self._parse_cpd_pred(cpd_pred)
        self._column_list = ['p', 'p_id', 'p_name', 'p_equation', 'p_genes']
        gene_map = self._reversible_map(gene_map)
        self._reaction_map_p = map_model_reactions(
            self._model1, self._model2, self._cpd_map, self._cpd_score, nproc,
            outpath, log=log, gene=gene, compartment_map=compartment_map,
            gene_map=gene_map)

    @property
    def model1(self):
        return self._model1

    @property
    def model2(self):
        return self._model2

    def map(self, r1, r2):
        return self._reaction_map_p[0][r1, r2]

    def _reversible_map(self, genemap):
        newmap = {}
        for k, v in genemap.items():
            newmap[k] = v
            newmap[v] = k
        return newmap

    def _parse_cpd_pred(self, cpd_pred):
        self._cpd_map = defaultdict(set)
        self._cpd_score = dict()
        for pair, score in cpd_pred.items():
            self._cpd_map[pair[0]].add(pair[1])
            self._cpd_score[pair[0]] = score

    def get_raw_map(self):
        """Return pandas.DataFrame style of raw mapping table."""
        reaction_result = pd.DataFrame({
            self._column_list[i]: self._reaction_map_p[i]
            for i in range(len(self._column_list))
        })
        return reaction_result

    def get_best_map(self, threshold_reaction=0):
        """Return pandas.DataFrame style of best mapping for each query."""
        raw_map = self.get_raw_map()
        query = [i for i, j in raw_map.index.values]
        # use rank instead of idxmax to output multiple top-hitts
        best_index = raw_map.iloc[:, 0].groupby(query).rank(
            method='min', ascending=False).loc[lambda x: x == 1].index.values
        reaction_best = raw_map.loc[best_index]
        reaction_best.query(
            'p >= @threshold_reaction', inplace=True)
        return reaction_best


def compound_id_likelihood(c1, c2, compound_prior, compound_id_marg):
    if util.id_equals(c1.id, c2.id):
        p_match = 0.4
        p_marg = compound_id_marg
        p_no_match = max(
            0, (p_marg - p_match * compound_prior) / (1.0 - compound_prior))
    else:
        p_match = 0.6
        p_marg = 1.0 - compound_id_marg
        p_no_match = max(
            0, (p_marg - p_match * compound_prior) / (1.0 - compound_prior))

    return p_match, p_no_match


def compound_name_likelihood(c1, c2, compound_prior, compound_name_marg):
    if util.name_equals(c1.name, c2.name):
        p_match = 0.60
        p_marg = compound_name_marg
        p_no_match = max(
            0, (p_marg - p_match * compound_prior) / (1.0 - compound_prior))
    else:
        p_match = 0.40
        p_marg = 1.0 - compound_name_marg
        p_no_match = max(
            0, (p_marg - p_match * compound_prior) / (1.0 - compound_prior))

    return p_match, p_no_match  # , p_marg


def compound_charge_likelihood(
        c1, c2, compound_prior,
        compound_charge_equal_marg, compound_charge_not_equal_marg):
    if c1.charge is None or c2.charge is None:
        # p value of observing undefined charge
        # it is independent of the condition of match or not
        p_match = 1
        p_no_match = 1
    elif c1.charge == c2.charge:
        p_match = 0.9
        p_no_match = max(
            0,
            ((compound_charge_equal_marg - p_match * compound_prior) /
             (1.0 - compound_prior)))
    else:
        p_match = 0.1
        p_no_match = max(
            0,
            ((compound_charge_not_equal_marg - p_match * compound_prior) /
             (1.0 - compound_prior)))

    return p_match, p_no_match


def compound_formula_likelihood(
        c1, c2, compound_prior,
        compound_formula_equal_marg, compound_formula_not_equal_marg):
    if c1.formula is None or c2.formula is None:
        # p value of observing undefined formula
        # it is independent of the condition of match or not
        p_match = 1
        p_no_match = 1
    elif util.formula_equals(c1.formula, c2.formula, c1.charge, c2.charge):
        p_match = 0.9
        p_no_match = max(
            0,
            ((compound_formula_equal_marg - p_match * compound_prior) /
             (1.0 - compound_prior)))
    else:
        p_match = 0.1
        p_no_match = max(
            0,
            ((compound_formula_not_equal_marg - p_match * compound_prior) /
             (1.0 - compound_prior)))

    return p_match, p_no_match


def compound_kegg_likelihood(
        c1, c2, compound_prior,
        compound_kegg_equal_marg, compound_kegg_not_equal_marg):
    if c1.kegg is None or c2.kegg is None:
        # p value of observing undefined KEGG
        # it is independent of the condition of match or not
        p_match = 1
        p_no_match = 1
    elif c1.kegg == c2.kegg:
        p_match = 0.65
        p_no_match = max(
            0,
            ((compound_kegg_equal_marg - p_match * compound_prior) /
             (1.0 - compound_prior)))
    else:
        p_match = 0.35
        p_no_match = max(
            0,
            ((compound_kegg_not_equal_marg - p_match * compound_prior) /
             (1.0 - compound_prior)))

    return p_match, p_no_match


def reaction_id_likelihood(
        r1, r2, reaction_prior,
        reaction_id_equal_marg, reaction_id_not_equal_marg):
    if util.id_equals(r1.id, r2.id):
        p_match = 0.3
        p_no_match = max(
            0,
            ((reaction_id_equal_marg - p_match * reaction_prior) /
             (1.0 - reaction_prior)))
    else:
        p_match = 0.7
        p_no_match = max(
            0,
            ((reaction_id_not_equal_marg - p_match * reaction_prior) /
             (1.0 - reaction_prior)))

    return p_match, p_no_match


def reaction_name_likelihood(r1, r2, reaction_prior, reaction_name_marg):
    if util.name_equals(r1.name, r2.name):
        p_match = 0.2
        p_no_match = max(
            0,
            ((reaction_name_marg - p_match * reaction_prior) /
             (1.0 - reaction_prior)))
    else:
        p_match = 0.8
        p_no_match = max(
            0,
            ((1.0 - reaction_name_marg - p_match * reaction_prior) /
             (1.0 - reaction_prior)))

    return p_match, p_no_match


def reaction_equation_mapping_approx_max_likelihood(
        cpd_set1, cpd_set2, cpd_map, cpd_score, compartment_map={}):
    """Calculate equation likelihood based on compound mapping."""
    p_match = 1.0
    p_no_match = 1.0

    compartment_dict = defaultdict(set)
    for c in cpd_set2:
        compartment_dict[c.name].add(c.compartment)

    # get the possible best-match pairs, score as the key
    best_match = defaultdict(set)
    for c1 in cpd_set1:
        for c2 in cpd_map[c1.name].intersection(compartment_dict.keys()):
            # the compound pair should be in the same compartment
            if (compartment_map.get(c1.compartment, c1.compartment)
                    in compartment_dict[c2]):
                best_match[cpd_score[c1.name]].add((c1.name, c2))

    # translate each element to compound id
    cpd_set1 = set(c.name for c in cpd_set1)
    cpd_set2 = set(c.name for c in cpd_set2)

    # sort by scores, get the most possible compound pairs first
    for score in sorted(best_match.keys(), reverse=True):
        for c1, c2 in best_match[score]:
            if (c1 in cpd_set1) and (c2 in cpd_set2):
                # the possibility that compounds are equal
                p_match *= score * 0.9 + (1 - score) * 0.1
                p_no_match *= score * 0.1 + (1 - score) * 0.9
                cpd_set1.remove(c1)
                cpd_set2.remove(c2)
                if len(cpd_set1) == 0 or len(cpd_set2) == 0:
                    break

    for c in cpd_set1:
        p_match *= 0.1
        p_no_match *= 0.9
    for c in cpd_set2:
        p_match *= 0.1
        p_no_match *= 0.9

    return p_match, p_no_match


def reaction_equation_compound_mapping_likelihood(
        r1, r2, *args, **kwargs):
    """Get the likelihood of reaction equations

    Args:
        r1, r2: two `RactionEntry` objects to be compared
    args, kwargs:
        cpd_map: dictionary mapping compound id in the query model to a set
                 of best-mapping compound ids in the target model.
        cpd_score: dictionary mapping compound id in the query model to
                   its best mapping score during compound mapping.
        compartment_map: dictionary mapping compartment id in the query model
                         to the id in the target model.
    """
    if r1.equation is None or r2.equation is None:
        # p value of observing undefined equation
        # it is independent of the condition of match or not
        p_match = 1
        p_no_match = 1
    else:
        p_match, p_no_match = get_best_p_value_set(r1, r2, *args, **kwargs)

    return p_match, p_no_match


def get_best_p_value_set(r1, r2, *args, **kwargs):
    """Assume equations may have reversed direction, report best mapping p."""
    cpd_set1_left = get_cpd_set(r1.equation, left=True)
    cpd_set1_right = get_cpd_set(r1.equation, left=False)
    cpd_set2_left = get_cpd_set(r2.equation, left=True)
    cpd_set2_right = get_cpd_set(r2.equation, left=False)

    # assume equations have the same direction
    p_forward_match, p_forward_no_match = merge_partial_p_set(
        cpd_set1_left, cpd_set2_left,
        cpd_set1_right, cpd_set2_right, *args, **kwargs)
    # assume equations have the reversed direction
    p_reverse_match, p_reverse_no_match = merge_partial_p_set(
        cpd_set1_left, cpd_set2_right,
        cpd_set1_right, cpd_set2_left, *args, **kwargs)

    # maintain the direction with better p values
    if (p_forward_match / p_forward_no_match >=
            p_reverse_match / p_reverse_no_match):
        p_match = p_forward_match
        p_no_match = p_forward_no_match
    else:
        p_match = p_reverse_match
        p_no_match = p_reverse_no_match

    return p_match, p_no_match


def merge_partial_p_set(cpd_set1_left, cpd_set2_left,
                        cpd_set1_right, cpd_set2_right, *args, **kwargs):
    """Merge the left hand side and right hand side p values together.

    The compound mapping is done separately on left hand side and
    right hand side.
    Then the corresponding p_match and p_no_match are merged together.
    """
    p_set_left = \
        reaction_equation_mapping_approx_max_likelihood(
            cpd_set1_left, cpd_set2_left, *args, **kwargs)
    p_set_right = \
        reaction_equation_mapping_approx_max_likelihood(
            cpd_set1_right, cpd_set2_right, *args, **kwargs)
    p_match = p_set_left[0] * p_set_right[0]
    p_no_match = p_set_left[1] * p_set_right[1]
    return p_match, p_no_match


def get_cpd_set(equation, left=True):
    if left:  # value of left-side compound is negtive
        coef = 1
    else:  # value of right-side compound is positive
        coef = -1
    # pick compounds at one side only
    cpd_set = set(
        compound
        for compound, value in equation.compounds
        if coef * value < 0)
    return cpd_set


def reaction_genes_likelihood(r1, r2, reaction_prior, reaction_genes_marg,
                              reaction_genes_not_equal_marg, gene_map={}):
    if r1.genes is None or r2.genes is None:
        p_match = 1
        p_no_match = 1
    elif util.genes_equals(r1.genes, r2.genes, gene_map):
        p_match = 0.2
        p_no_match = max(
            0,
            ((reaction_genes_marg - p_match * reaction_prior) /
             (1.0 - reaction_prior)))
    else:
        p_match = 0.8
        p_no_match = max(
            0,
            ((reaction_genes_not_equal_marg - p_match * reaction_prior) /
             (1.0 - reaction_prior)))

    return p_match, p_no_match


def fake_likelihood(e1, e2):
    """Generate fake likelihood if corresponding mapping is not required."""
    return 1, 1


def generate_likelihood(tasks):
    pair, likelihood, args, kwargs = tasks
    e1, e2 = pair
    p1, p2 = likelihood(e1, e2, *args, **kwargs)
    return e1.id, e2.id, p1, p2


def pairwise_likelihood(pool, chunksize, model1, model2, likelihood,
                        *args, **kwargs):
    """Compute likelihood of all pairwise comparisons.

    Returns likelihoods as a dataframe with a column for each hypothesis.
    """
    tasks = (((e1, e2), likelihood, args, kwargs)
             for e1, e2 in product(itervalues(model1), itervalues(model2)))
    result = pool.map(generate_likelihood, tasks, chunksize=chunksize)
    return pd.DataFrame.from_records(result, index=('e1', 'e2'),
                                     columns=('e1', 'e2', 'p1', 'p2'))


def likelihood_products(likelihood_dfs):
    """Combine likelihood dataframes."""
    return reduce(operator.mul, likelihood_dfs, 1.0)


def bayes_posterior(prior, likelihood_df):
    """Calculate posterior given likelihoods and prior."""
    p_1 = prior * likelihood_df.iloc[:, 0]
    p_2 = (1.0 - prior) * likelihood_df.iloc[:, 1]
    return p_1 / (p_1 + p_2)


def parallel_equel(tasks):
    func, params = tasks
    return func(*params)


def map_model_compounds(model1, model2, nproc=1, outpath='.',
                        log=False, kegg=False):
    """Map compounds of two models."""
    compound_pairs = len(model1.compounds) * len(model2.compounds)

    # Compound prior
    # For the prior, use a guesstimate that 95% of the
    # smaller model can be mapped.
    compound_prior = (0.95 * min(len(model1.compounds),
                                 len(model2.compounds))) / compound_pairs

    # Initialize parallel pool of workers
    chunksize = compound_pairs // nproc
    pool = Pool(nproc)

    t = time.time()
    # Compound ID
    print('Calculating compound ID likelihoods...', end=' ')
    sys.stdout.flush()

    # Marginal probability of observing two equal compound IDs
    tasks = ((util.id_equals, (c1.id, c2.id)) for c1, c2 in product(
        itervalues(model1.compounds), itervalues(model2.compounds)))
    result = pool.map(parallel_equel, tasks, chunksize=chunksize)
    compound_id_marg = sum(result) / float(compound_pairs)

    compound_id_likelihoods = pairwise_likelihood(
        pool, chunksize, model1.compounds, model2.compounds,
        compound_id_likelihood, compound_prior, compound_id_marg)

    print('%.2f seconds' % (time.time() - t))
    t = time.time()
    # Compound name
    print('Calculating compound name likelihoods...', end=' ')
    sys.stdout.flush()

    # Marginal probability of observing two similar names
    tasks = ((util.name_equals, (c1.name, c2.name)) for c1, c2 in product(
        itervalues(model1.compounds), itervalues(model2.compounds)))
    result = pool.map(parallel_equel, tasks, chunksize=chunksize)
    compound_name_marg = sum(result) / float(compound_pairs)

    compound_name_likelihoods = pairwise_likelihood(
        pool, chunksize, model1.compounds, model2.compounds,
        compound_name_likelihood, compound_prior, compound_name_marg)

    print('%.2f seconds' % (time.time() - t))
    t = time.time()

    # Compound charge
    print('Calculating compound charge likelihoods...', end=' ')
    sys.stdout.flush()

    # Marginal probability of observing two compounds with the same charge
    compound_charge_equal_marg = sum(
        c1.charge is not None and
        c2.charge is not None and
        c1.charge == c2.charge
        for c1, c2 in product(
            itervalues(model1.compounds), itervalues(model2.compounds))
    ) / compound_pairs

    # Marginal probability of observing two compounds with different charge
    compound_charge_not_equal_marg = sum(
        c1.charge is not None and
        c2.charge is not None and
        c1.charge != c2.charge
        for c1, c2 in product(
            itervalues(model1.compounds), itervalues(model2.compounds))
    ) / compound_pairs

    compound_charge_likelihoods = pairwise_likelihood(
        pool, chunksize, model1.compounds, model2.compounds,
        compound_charge_likelihood,
        compound_prior,
        compound_charge_equal_marg,
        compound_charge_not_equal_marg)

    print('%.2f seconds' % (time.time() - t))
    t = time.time()

    # Compound formula
    print('Calculating compound formula likelihoods...', end=' ')
    sys.stdout.flush()

    # Marginal probability of observing two compounds with the same formula
    tasks = ((
        util.formula_equals,
        (c1.formula, c2.formula, c1.charge, c2.charge))
        for c1, c2 in product(
            itervalues(model1.compounds), itervalues(model2.compounds)))
    result = pool.map(parallel_equel, tasks, chunksize=chunksize)
    compound_formula_equal_marg = sum(result) / float(compound_pairs)

    # Marginal probability of observing two compounds with different formula
    compound_formula_not_equal_marg = 1.0 - compound_formula_equal_marg - (
        sum(c1.formula is None or c2.formula is None
            for c1, c2 in product(itervalues(model1.compounds),
                                  itervalues(model2.compounds))) /
        compound_pairs)

    compound_formula_likelihoods = pairwise_likelihood(
        pool, chunksize, model1.compounds, model2.compounds,
        compound_formula_likelihood,
        compound_prior, compound_formula_equal_marg,
        compound_formula_not_equal_marg)

    print('%.2f seconds' % (time.time() - t))
    t = time.time()

    # Compound KEGG id
    if kegg:  # run KEGG id mapping
        print('Calculating compound KEGG ID likelihoods...', end=' ')
        sys.stdout.flush()

        # Marginal probability of observing two compounds
        # where KEGG ids are equal
        compound_kegg_equal_marg = sum(
            c1.kegg is not None and
            c2.kegg is not None and
            c1.kegg == c2.kegg
            for c1, c2 in product(
                itervalues(model1.compounds),
                itervalues(model2.compounds))
        ) / compound_pairs

        # Marginal probability of observing two compounds
        # where KEGG ids are different
        compound_kegg_not_equal_marg = sum(
            c1.kegg is not None and
            c2.kegg is not None and
            c1.kegg != c2.kegg for c1, c2 in product(
                itervalues(model1.compounds),
                itervalues(model2.compounds))
        ) / compound_pairs

        compound_kegg_likelihoods = pairwise_likelihood(
            pool, chunksize, model1.compounds, model2.compounds,
            compound_kegg_likelihood,
            compound_prior, compound_kegg_equal_marg,
            compound_kegg_not_equal_marg)

        print('%.2f seconds' % (time.time() - t))
        t = time.time()
    else:  # run fake mapping
        compound_kegg_likelihoods = pairwise_likelihood(
            pool, chunksize, model1.compounds, model2.compounds,
            fake_likelihood)

    pool.close()
    pool.join()

    if log:
        merge_result = pd.merge(compound_id_likelihoods,
                                compound_name_likelihoods,
                                left_index=True, right_index=True,
                                suffixes=('_id', '_name'))
        merge_result = pd.merge(merge_result, compound_charge_likelihoods,
                                left_index=True, right_index=True,
                                suffixes=('_name', '_charge'))
        merge_result = pd.merge(merge_result, compound_formula_likelihoods,
                                left_index=True, right_index=True,
                                suffixes=('_charge', '_formula'))
        merge_result = pd.merge(merge_result, compound_kegg_likelihoods,
                                left_index=True, right_index=True,
                                suffixes=('_formula', '_kegg'))

        merge_result.to_csv(outpath + '/compound_log.tsv', sep='\t')

    all_likelihoods = [compound_id_likelihoods,
                       compound_name_likelihoods,
                       compound_charge_likelihoods,
                       compound_formula_likelihoods,
                       compound_kegg_likelihoods]

    return (bayes_posterior(compound_prior,
                            likelihood_products(all_likelihoods)),
            bayes_posterior(compound_prior, compound_id_likelihoods),
            bayes_posterior(compound_prior, compound_name_likelihoods),
            bayes_posterior(compound_prior, compound_charge_likelihoods),
            bayes_posterior(compound_prior, compound_formula_likelihoods),
            bayes_posterior(compound_prior, compound_kegg_likelihoods))


def map_model_reactions(model1, model2, cpd_map, cpd_score, nproc=1,
                        outpath='.', log=False, gene=False, compartment_map={},
                        gene_map={}):
    """Map reactions of two models."""
    # Mapping of reactions
    reaction_pairs = len(model1.reactions) * len(model2.reactions)

    # Reaction prior
    # For the prior, use a guesstimate that 95%
    # of the smaller model can be mapped.
    reaction_prior = (0.95 * min(len(model1.reactions),
                                 len(model2.reactions))) / reaction_pairs

    # Initialize parallel pool of workers
    chunksize = reaction_pairs // nproc
    pool = Pool(nproc)

    # Reaction ID
    t = time.time()
    print('Calculating reaction ID likelihoods...', end=' ')
    sys.stdout.flush()

    # Marginal probability of observing two reactions with the same ids.
    tasks = ((util.id_equals, (r1.id, r2.id)) for r1, r2 in product(
        itervalues(model1.reactions),
        itervalues(model2.reactions)))
    result = pool.map(parallel_equel, tasks, chunksize=chunksize)
    reaction_id_equal_marg = sum(result) / float(reaction_pairs)

    # Marginal probability of observing two reactions with different ids.
    reaction_id_not_equal_marg = 1.0 - reaction_id_equal_marg

    reaction_id_likelihoods = pairwise_likelihood(
        pool, chunksize, model1.reactions, model2.reactions,
        reaction_id_likelihood,
        reaction_prior, reaction_id_equal_marg, reaction_id_not_equal_marg)

    print('%.2f seconds' % (time.time() - t))
    t = time.time()
    # Reaction name
    print('Calculating reaction name likelihoods...', end=' ')
    sys.stdout.flush()

    # Marginal probability of observing two reactions with the same name.
    tasks = ((util.name_equals, (r1.name, r2.name)) for r1, r2 in product(
        itervalues(model1.reactions),
        itervalues(model2.reactions)))
    result = pool.map(parallel_equel, tasks, chunksize=chunksize)
    reaction_name_equal_marg = sum(result) / float(reaction_pairs)

    reaction_name_likelihoods = pairwise_likelihood(
        pool, chunksize, model1.reactions, model2.reactions,
        reaction_name_likelihood, reaction_prior, reaction_name_equal_marg)

    print('%.2f seconds' % (time.time() - t))
    t = time.time()

    # Reaction equation

    print('Calculating reaction equation likelihoods...', end=' ')
    sys.stdout.flush()
    reaction_equation_likelihoods = pairwise_likelihood(
        pool, chunksize, model1.reactions, model2.reactions,
        reaction_equation_compound_mapping_likelihood,
        cpd_map, cpd_score, compartment_map)

    print('%.2f seconds' % (time.time() - t))
    t = time.time()

    # Reaction genes
    # For each gene, the marginal probability of observing that gene
    # in each model. We use this as an approximation of the probability of
    # observing a pair of genes in two reactions given that the reaction
    # do _not_ match.
    if gene:
        print('Calculating reaction genes likelihoods...', end=' ')
        sys.stdout.flush()

        # Marginal probability of observing two reactions with
        # equal gene associations.
        tasks = ((util.genes_equals, (r1.genes, r2.genes))
                 for r1, r2 in product(
                     itervalues(model1.reactions),
                     itervalues(model2.reactions)))
        result = pool.map(parallel_equel, tasks, chunksize=chunksize)
        reaction_genes_equal_marg = sum(result) / float(reaction_pairs)

        # Marginal probability of observing two reactions with unequal
        # gene associations.
        reaction_genes_not_equal_marg = 1.0 - reaction_genes_equal_marg - (
            sum(r1.genes is None or r2.genes is None
                for r1, r2 in product(itervalues(model1.reactions),
                                      itervalues(model2.reactions))) /
            reaction_pairs)

        reaction_genes_likelihoods = pairwise_likelihood(
            pool, chunksize, model1.reactions, model2.reactions,
            reaction_genes_likelihood,
            reaction_prior, reaction_genes_equal_marg,
            reaction_genes_not_equal_marg, gene_map)

        print('%.2f seconds' % (time.time() - t))
        t = time.time()
    else:
        reaction_genes_likelihoods = pairwise_likelihood(
            pool, chunksize, model1.reactions, model2.reactions,
            fake_likelihood)

    pool.close()
    pool.join()

    if log:
        merge_result = pd.merge(reaction_id_likelihoods,
                                reaction_name_likelihoods,
                                left_index=True, right_index=True,
                                suffixes=('_id', '_name'))
        merge_result = pd.merge(merge_result, reaction_equation_likelihoods,
                                left_index=True, right_index=True,
                                suffixes=('_name', '_equation'))
        merge_result = pd.merge(merge_result, reaction_genes_likelihoods,
                                left_index=True, right_index=True,
                                suffixes=('_equation', '_genes'))

        merge_result.to_csv(outpath + '/reaction_log.tsv', sep='\t')

    all_likelihoods = [reaction_id_likelihoods,
                       reaction_name_likelihoods,
                       reaction_equation_likelihoods,
                       reaction_genes_likelihoods]

    return (bayes_posterior(reaction_prior,
                            likelihood_products(all_likelihoods)),
            bayes_posterior(reaction_prior, reaction_id_likelihoods),
            bayes_posterior(reaction_prior, reaction_name_likelihoods),
            bayes_posterior(reaction_prior, reaction_equation_likelihoods),
            bayes_posterior(reaction_prior, reaction_genes_likelihoods))


def check_cpd_charge(compound, source):
    if compound.charge is not None:
        if isinstance(compound.charge, int):
            return True
        else:
            logger.warning(
                "Compound charge should be an integer, however, charge of "
                "compound '{}' in '{}' model is '{}', which is invalid. "
                "Please remove or fix it before running modelmapping "
                "command.".format(compound.id, source, compound.charge))
            return False
