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

# -*- coding: utf-8 -*-

"""Generate a common model for two distinct metabolic models."""

from __future__ import print_function, division
from future.utils import itervalues

import re
from collections import namedtuple
from builtins import object, input
import operator
from itertools import product
import time
import sys

from psamm.formula import Formula
import psamm.bayesian as bayesian
import psamm.translate_id as tr_id
import psamm.manual_curation as curation
from psamm.util import mkdir_p
from psamm.datasource.native import ModelReader
from psamm.command import Command


CompoundEntry = namedtuple(
    'CompoundEntry',
    ['id', 'name', 'formula', 'charge', 'kegg', 'cas'])

ReactionEntry = namedtuple(
    'ReactionEntry',
    ['id', 'name', 'genes', 'equation', 'subsystem', 'ec'])


class MetabolicModel(object):
    def __init__(self, name, compounds, reactions):
        self._name = name
        self._read_compounds(compounds)
        self._read_reactions(reactions)

        self._genes = set()
        for r in itervalues(self._reactions):
            if r.genes is not None:
                self._genes.update(r.genes)

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
            if genes is not None:
                # Delete operators in string
                genes = re.sub(r'[,\(\)(or)(and)]*', '', genes)
                genes = re.split(r'\s+', genes)  # Transfer string to list
                genes = frozenset(genes)
            else:
                genes = None

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

    def check_reaction_compounds(self):
        """Check that reaction compounds are defined in the model"""
        print('Checking model: {}'.format(self.name))
        consistent = True
        for reaction in itervalues(self.reactions):
            if reaction.equation is not None:
                for compound, value in reaction.equation.compounds:
                    if compound.name not in self.compounds:
                        print((
                            '{} in reaction {} '
                            'is not in compound list').format(
                                compound.name, reaction.id))
                        consistent = False
        return consistent


def read_mapping_file(f):
    """Read a mapping file and yield source, target-set tuples"""
    for line in f:
        source, target = line.strip().split(None, 1)
        if target == '-':
            target_set = None
        else:
            target_set = set(s for s in target.split(','))
        yield source, target_set


class Confusion(object):
    """Confusion matrix."""

    def __init__(self, tp=0, tn=0, fp=0, fn=0):
        self.tp, self.tn, self.fp, self.fn = tp, tn, fp, fn

    @property
    def positive(self):
        return self.tp + self.fn

    @property
    def negative(self):
        return self.fp + self.tn

    @property
    def total(self):
        return self.tp + self.tn + self.fp + self.fn

    @property
    def tp_rate(self):
        return self.tp if self.tp == 0 else self.tp / self.positive

    @property
    def fp_rate(self):
        return self.fp if self.fp == 0 else self.fp / self.negative

    @property
    def specificity(self):
        return self.tn if self.tn == 0 else self.tn / self.negative

    @property
    def precision(self):
        return self.tp if self.tp == 0 else self.tp / (self.tp + self.fp)

    @property
    def accuracy(self):
        return (self.tp + self.tn) / self.total

    @property
    def balanced_accuracy(self):
        return (self.tp_rate + self.specificity) / 2.0

    @property
    def informedness(self):
        return self.tp_rate + self.specificity - 1.0

    def __repr__(self):
        return 'Confusion(tp={}, tn={}, fp={}, fn={})'.format(
            self.tp, self.tn, self.fp, self.fn)


def evaluate_roc_points(model1, model2, predictor, actual):
    """Return confusion matrix for each point on the ROC curve"""

    # Count total positive and negative pairs
    positive, negative = 0, 0
    for e1, e2 in product(model1, model2):
        correct = e1 in actual and actual[e1] is not None and e2 in actual[e1]
        positive += correct
        negative += not correct

    tp, fp = 0, 0
    p_prev = None
    for pair, p in sorted((((e1, e2), predictor.map(e1, e2))
                           for e1, e2 in product(model1, model2)),
                          key=operator.itemgetter(1), reverse=True):
        e1, e2 = pair
        if p != p_prev:
            yield p, Confusion(
                tp=tp, fp=fp, fn=positive - tp, tn=negative - fp)
            p_prev = p
        correct = e1 in actual and actual[e1] is not None and e2 in actual[e1]
        tp += correct
        fp += not correct

    cf = Confusion(tp=tp, fp=fp, fn=positive - tp, tn=negative - fp)
    yield -float('inf'), cf


def write_roc_curve(out_file, model1, model2, predictor, actual):
    """Calculate ROC points and write them to file"""
    for p, conf in evaluate_roc_points(model1, model2,
                                       predictor, actual):
        out_file.write('{}\t{}\t{}\t{}\t{}\n'.format(
            p, conf.tp, conf.tn, conf.fp, conf.fn))


class ModelMappingCommand(Command):
    """Generate a common model for two distinct metabolic models."""

    @classmethod
    def init_parser(cls, parser):
        """Initialize argument parser"""
        subparsers = parser.add_subparsers(title="Tools")

        # model mapping subcommand
        parser_mm = subparsers.add_parser(
            'map', help='Bayesian model mapping')
        parser_mm.set_defaults(which='mm')
        parser_mm.add_argument(
            '--dest-model', type=str, required=True, help='Model to map to')
        parser_mm.add_argument(
            '--compound-map', type=str, nargs='?',
            help='Actual compound mapping, evaluation only')
        parser_mm.add_argument(
            '--reaction-map', type=str, nargs='?',
            help='Actual reaction mapping, evaluation only')
        parser_mm.add_argument(
            '-c', '--consistency-check', action='store_true',
            help=(
                'Check the consistency of compound id '
                'between compounds.yaml and reactions.yaml'))
        parser_mm.add_argument(
            '-n', '--nproc', type=int, action='store', default=1,
            help='Number of processes to use (default = 1)')
        parser_mm.add_argument(
            '-o', '--outpath', type=str,
            action='store', default='.',
            help=(
                'Path to store output files '
                '(default = current position)'))
        parser_mm.add_argument(
            '--threshold-compound', type=float, default=0,
            help=(
                'likelihood threshold for compound mapping '
                '(default = 0)'))
        parser_mm.add_argument(
            '--threshold-reaction', type=float, default=0,
            help=(
                'likelihood threshold for reaction mapping '
                '(default = 0)'))
        parser_mm.add_argument(
            '--map-compound-kegg', action='store_true',
            help='Compare KEGG id during compound mapping')
        parser_mm.add_argument(
            '--map-reaction-gene', action='store_true',
            help='Compare gene list during reaction mapping')
        parser_mm.add_argument(
            '--raw', action='store_true',
            help=(
                'Output raw pairwise mapping results as well, '
                'may be very big with large models'))
        parser_mm.add_argument(
            '--log', action='store_true',
            help=(
                'Output log file of the p_match and p_no_match '
                'for each feature, '
                'may be very big with large models'))

        # manual curation subcommand
        parser_c = subparsers.add_parser(
            'manual_curation',
            help='Interactive mode to check mapping result')
        parser_c.set_defaults(which='curation')
        parser_c.add_argument(
            '--dest-model', type=str, required=True, help='Model to map to')
        parser_c.add_argument(
            '--compound-map', type=str, required=True,
            help='Compound model mapping result')
        parser_c.add_argument(
            '--reaction-map', type=str, required=True,
            help='Reaction model mapping result')
        parser_c.add_argument(
            '--curated-compound-map', type=str, required=True,
            help=(
                'File to store curated compound mapping, if file exists, '
                'resume previous curation'))
        parser_c.add_argument(
            '--curated-reaction-map', type=str, required=True,
            help=(
                'File to store curated reaction mapping, if file exists, '
                'resume previous curation'))

        # translate id subcommand
        parser_t = subparsers.add_parser(
            'translate_id', help='Translate ids based on mapping result')
        parser_t.set_defaults(which='translateid')
        parser_t.add_argument(
            '--compound-map', type=str, required=True,
            help=(
                'Tab-delimited table, the first two columns store the '
                'original and target ids'))
        parser_t.add_argument(
            '--reaction-map', type=str, required=True,
            help=(
                'Tab-delimited table, the first two columns store the '
                'original and target ids'))
        parser_t.add_argument(
            '-o', '--outpath', type=str,
            action='store', default='.',
            help=(
                'Path to store output files '
                '(default = current position)'))

    def run(self):
        """Run model mapping command."""
        which_command = self._args.which
        if which_command == 'mm':
            self._model_mapping()
        if which_command == 'translateid':
            self._translate_id()
        if which_command == 'curation':
            self._curation()

    def _model_mapping(self):
        """Run model mapping"""
        # Parse models
        model1 = ModelReader.reader_from_path(self._args.model)
        model2 = ModelReader.reader_from_path(self._args.dest_model)

        # Read model into internal format
        model1 = MetabolicModel(
            model1.name, model1.parse_compounds(), model1.parse_reactions())
        model2 = MetabolicModel(
            model2.name, model2.parse_compounds(), model2.parse_reactions())

        # Model summaries
        model1.print_summary()
        print('\n')
        model2.print_summary()
        print('\n')

        # Load actual model mappings
        actual_compound_mapping = None
        if self._args.compound_map is not None:
            with open(self._args.compound_map, 'r') as f:
                actual_compound_mapping = dict(read_mapping_file(f))

        actual_reaction_mapping = None
        if self._args.reaction_map is not None:
            with open(self._args.reaction_map, 'r') as f:
                actual_reaction_mapping = dict(read_mapping_file(f))

        # Check models
        if (
            not (
                model1.check_reaction_compounds() and
                model2.check_reaction_compounds()) and
                self._args.consistency_check):
            quit((
                '\nError: '
                'equations have something not listed in compounds.yaml, '
                'please check it!'))

        mkdir_p(self._args.outpath)

        if (actual_compound_mapping is not None):
            mkdir_p(self._args.outpath + '/roc')

        # Bayesian classifier
        print('Using %i processes...' % (self._args.nproc))
        t = time.time()
        cpd_bayes_pred = bayesian.BayesianCompoundPredictor(
            model1, model2, self._args.nproc, self._args.outpath,
            log=self._args.log, kegg=self._args.map_compound_kegg)
        print(
            'It took %s seconds to calculate compound mapping...'
            % (time.time() - t))

        print('Writing output...')
        sys.stdout.flush()
        t = time.time()
        # Write out ROC curve results
        if (actual_compound_mapping is not None):
            with open(self._args.outpath +
                      '/roc/compound_bayes.tsv', 'w') as f:
                write_roc_curve(f, cpd_bayes_pred.model1.compounds,
                                cpd_bayes_pred.model2.compounds,
                                cpd_bayes_pred,
                                actual_compound_mapping)

        # Write out classifier results

        # if (actual_compound_mapping is not None):
        #     for c1, c2 in product(cpd_bayes_pred.model1.compounds,
        #                           cpd_bayes_pred.model2.compounds):
        #         maps = (c1 in actual_compound_mapping and
        #                 actual_compound_mapping[c1] is not None and
        #                 c2 in actual_compound_mapping[c1])
        #         p = cpd_bayes_pred.map(c1, c2)
        #         if (p > self._args.threshold_compound and not maps):
        #             print('Incorrect map: {}\t{}\t{}'.format(c2, c1, p))
        #         if (p <= self._args.threshold_compound and maps):
        #             print('Incorrect unmap: {}\t{}\t{}'.format(c2, c1, p))

        # Parse and output raw mapping
        if self._args.raw:
            cpd_bayes_pred.get_raw_map().to_csv(
                self._args.outpath + '/bayes_compounds.tsv', sep='\t')

        # Output best mapping
        compound_best = cpd_bayes_pred.get_best_map(
            self._args.threshold_compound)
        compound_best.to_csv(
            self._args.outpath + '/bayes_compounds_best.tsv', sep='\t')
        print('It took %s seconds to write output...' % (time.time() - t))
        sys.stdout.flush()

        # Bayesian classifier
        t = time.time()
        rxn_bayes_pred = bayesian.BayesianReactionPredictor(
            model1, model2, compound_best.loc[:, 'p'],
            self._args.nproc, self._args.outpath,
            log=self._args.log, gene=self._args.map_reaction_gene)
        print(
            'It took %s seconds to calculate reaction mapping...'
            % (time.time() - t))

        print('Writing output...')
        sys.stdout.flush()
        t = time.time()
        # Write out ROC curve results
        if (actual_reaction_mapping is not None):
            with open(self._args.outpath +
                      '/roc/reaction_bayes.tsv', 'w') as f:
                write_roc_curve(f, rxn_bayes_pred.model1.reactions,
                                rxn_bayes_pred.model2.reactions,
                                rxn_bayes_pred,
                                actual_reaction_mapping)

        # Write out classifier results

        # if (actual_reaction_mapping is not None):
        #     for r1, r2 in product(rxn_bayes_pred.model1.reactions,
        #                         rxn_bayes_pred.model2.reactions):
        #         maps = (r1 in actual_reaction_mapping and
        #                 actual_reaction_mapping[r1] is not None and
        #                 r2 in actual_reaction_mapping[r1])
        #         p = rxn_bayes_pred.map(r1, r2)
        #         if (p > self._args.threshold_reaction and not maps):
        #             print('Incorrect map: {}\t{}\t{}'.format(r2, r1, p))
        #         if (p <= self._args.threshold_reaction and maps):
        #             print('Incorrect unmap: {}\t{}\t{}'.format(r2, r1, p))

        # Parse and output raw mapping
        if self._args.raw:
            rxn_bayes_pred.get_raw_map().to_csv(
                self._args.outpath + '/bayes_reactions.tsv', sep='\t')

        # Output best mapping
        reaction_best = rxn_bayes_pred.get_best_map(
            self._args.threshold_reaction)
        reaction_best.to_csv(
            self._args.outpath + '/bayes_reactions_best.tsv', sep='\t')
        print('It took %s seconds to write output...' % (time.time() - t))
        sys.stdout.flush()

    def _curation(self):
        # read models
        model1 = ModelReader.reader_from_path(self._args.model)
        model1 = model1.create_model()
        model2 = ModelReader.reader_from_path(self._args.dest_model)
        model2 = model2.create_model()

        # initiate curator
        curator = curation.Curator(
            self._args.compound_map, self._args.reaction_map,
            self._args.curated_compound_map, self._args.curated_reaction_map)

        print('Starting to do manual curation...\n')
        print('Type "stop" to stop\n')
        # iterate through reaction mapping
        for rmap in curator.reaction_map.iterrows():
            # skip already curated reactions
            if (curator.reaction_checked(rmap[0])
                    or curator.reaction_checked(rmap[0][0])):
                continue
            # check the compound mapping in current reaction
            compounds = curation.search_reaction(model1, [rmap[0][0]])
            dest_compounds = curation.search_reaction(model2, [rmap[0][1]])
            print('Below are the compound mapping involved:\n')
            for compound in compounds:
                for cmap in curator.compound_map.loc[compound].iterrows():
                    if (cmap is None
                            or cmap[0] not in dest_compounds
                            or curator.compound_checked((compound, cmap[0]))
                            or curator.compound_checked(compound)):
                        continue
                    print(cmap[1])
                    print('\n')
                    curation.search_compound(model1, [compound])
                    curation.search_compound(model2, [cmap[0]])
                    # waiting for legal curation input
                    while True:
                        ask = input(
                            ('True compound match? (y/n/ignore/save/stop, '
                             'type ignore to ignore this compound in future, '
                             'type save to save current progress, '
                             'type stop to save and exit): '))
                        if ask.lower() == 'n':
                            curator.add_mapping(
                                (compound, cmap[0]),
                                'c',
                                False
                            )
                            break
                        if ask.lower() == 'y':
                            curator.add_mapping(
                                (compound, cmap[0]),
                                'c',
                                True
                            )
                            break
                        if ask.lower() == 'ignore':
                            curator.add_ignore(compound, 'c')
                            break
                        if ask.lower() == 'stop':
                            curator.save()
                            exit(0)
                        if ask.lower() == 'save':
                            curator.save()

                    print('\n')
            print('Here is the reaction mapping:')
            print(rmap[1])
            print('\n')
            curation.search_reaction(model1, [rmap[0][0]])
            curation.search_reaction(model2, [rmap[0][1]])
            print(('These two reactions have the following curated compound '
                   'pairs:\n'))
            count = 0
            for compound in compounds:
                if compound in curator.curated_compound_map.index:
                    for cmap in curator.curated_compound_map.loc[
                            compound].iterrows():
                        if cmap[0] in dest_compounds:
                            print(compound, cmap[0], cmap[1]['p'])
                            count += 1
            print('\n')
            print('%i compounds in %s, ' % (len(compounds), rmap[0][0]))
            print('%i compounds in %s\n' % (len(dest_compounds), rmap[0][1]))
            print('%i curated compound pairs\n' % count)
            ask = ''
            # waiting for legal curation input
            while ask not in ['y', 'n', 'stop']:
                ask = input(
                    ('True reaction match? (y/n/save/stop, '
                     'type ignore to ignore this reaction in future, '
                     'type save to save current progress, '
                     'type stop to save and exit): ')
                )
                if ask.lower() == 'n':
                    curator.add_mapping(rmap[0], 'r', False)
                    break
                if ask.lower() == 'y':
                    curator.add_mapping(rmap[0], 'r', True)
                    break
                if ask.lower() == 'ignore':
                    curator.add_ignore(rmap[0][0], 'r')
                    break
                if ask.lower() == 'stop':
                    curator.save()
                    exit(0)
                if ask.lower() == 'save':
                    curator.save()
            print('\n', '=' * 50, '\n')
        curator.save()

    def _translate_id(self):
        """Translate ids"""
        # make output dir
        mkdir_p(self._args.outpath)

        # read mapping files
        cpd_mapping_id = tr_id.read_mapping(self._args.compound_map)
        rxn_mapping_id = tr_id.read_mapping(self._args.reaction_map)

        # make output model
        out_nm = tr_id.TranslatedModel(
            cpd_mapping_id,
            rxn_mapping_id,
            self._model)

        out_nm.write_model(dest=self._args.outpath, split_subsystem=False)
