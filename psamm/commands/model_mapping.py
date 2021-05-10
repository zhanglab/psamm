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
# Copyright 2020-2021  Elysha Sameth <esameth1@uri.edu>

# -*- coding: utf-8 -*-

"""Generate a common model for two distinct metabolic models."""

from __future__ import print_function, division

from builtins import object, input
import operator
from itertools import product
import time
import sys

import psamm.bayesian as bayesian
import psamm.translate_id as tr_id
import psamm.manual_curation as curation
from psamm.util import mkdir_p
from psamm.datasource.native import ModelReader
from psamm.command import Command
from psamm.datasource.reaction import Compound


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
        parser_mm.add_argument(
            '--compartment-map-file', type=str,
            help=(
                'A tab delimited file (.tsv) listing the compartment ids of '
                'query model in the first column, and the corresponding '
                'compartment ids of dest model in the second column. Required '
                'if the compartment ids in the two models are not consistent.'
            )
        )
        parser_mm.add_argument(
            '--gene-map-file', type=str,
            help=(
                'A tab delimited file (.tsv) listing the gene ids of '
                'query model in the first column, and the corresponding '
                'gene ids of dest model in the second column. Required '
                'if the gene ids in the two models are not consistent.'
            )
        )

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
        parser_c.add_argument(
            '--compartment-map', type=str,
            help=(
                'A tab delimited file (.tsv) listing the compartment ids of '
                'query model in the first column, and the corresponding '
                'compartment ids of dest model in the second column. Required '
                'if the compartment ids in the two models are not consistent.'
            )
        )
        parser_c.add_argument(
            '--compound-only', action='store_true',
            help='Set to only curate compound maps.')

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
            '--compartment-map', type=str, required=False,
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
        model1 = ModelReader.reader_from_path(
            self._args.model).create_model()
        model2 = ModelReader.reader_from_path(
            self._args.dest_model).create_model()

        # Check if any compound charge is not integer
        invalid_cpd = []
        for compound in model1.compounds:
            result = bayesian.check_cpd_charge(compound, model1.name)
            if not result:
                invalid_cpd.append(compound.id)
        for compound in model2.compounds:
            result = bayesian.check_cpd_charge(compound, model2.name)
            if not result:
                invalid_cpd.append(compound.id)
        if len(invalid_cpd) != 0:
            quit()

        # Read model into internal format
        model1 = bayesian.MappingModel(model1)
        model2 = bayesian.MappingModel(model2)

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

        # Load compartment mapping
        compartment_map = {}
        if self._args.compartment_map_file is not None:
            with open(self._args.compartment_map_file, 'r') as f:
                for row in f:
                    source, target = row.strip().split()
                    compartment_map[source] = target

        # Load gene mapping
        gene_map = {}
        if self._args.gene_map_file is not None:
            with open(self._args.gene_map_file, 'r') as f:
                for row in f:
                    source, target = row.strip().split()
                    gene_map[source] = target

        # Check models
        if (not (model1.check_reaction_compounds() and
                 model2.check_reaction_compounds())):
            quit((
                '\nError: '
                'equations have something not listed in compounds.yaml, '
                'please check it!'))

        mkdir_p(self._args.outpath)

        if actual_compound_mapping is not None:
            mkdir_p(self._args.outpath + '/roc')

        # Bayesian classifier
        print('Using %i processes...' % (self._args.nproc))
        t = time.time()
        cpd_bayes_pred = bayesian.BayesianCompoundPredictor(
            model1, model2, self._args.nproc, self._args.outpath,
            log=self._args.log, kegg=self._args.map_compound_kegg)
        print(
            'It took %.2f seconds to calculate compound mapping...'
            % (time.time() - t))

        print('Writing output...')
        sys.stdout.flush()
        t = time.time()
        # Write out ROC curve results
        if actual_compound_mapping is not None:
            with open(self._args.outpath +
                      '/roc/compound_bayes.tsv', 'w') as f:
                write_roc_curve(f, cpd_bayes_pred.model1.compounds,
                                cpd_bayes_pred.model2.compounds,
                                cpd_bayes_pred,
                                actual_compound_mapping)

        # Parse and output raw mapping
        if self._args.raw:
            cpd_bayes_pred.get_raw_map().to_csv(
                self._args.outpath + '/bayes_compounds.tsv', sep='\t',
                encoding='utf-8')

        # Output best mapping
        compound_best = cpd_bayes_pred.get_best_map(
            self._args.threshold_compound)
        compound_best.to_csv(
            self._args.outpath + '/bayes_compounds_best.tsv', sep='\t',
            encoding='utf-8')
        print('It took %.2f seconds to write output...' % (time.time() - t))
        sys.stdout.flush()

        # Bayesian classifier
        t = time.time()
        rxn_bayes_pred = bayesian.BayesianReactionPredictor(
            model1, model2, compound_best.loc[:, 'p'],
            self._args.nproc, self._args.outpath,
            log=self._args.log, gene=self._args.map_reaction_gene,
            compartment_map=compartment_map, gene_map=gene_map)
        print(
            'It took %.2f seconds to calculate reaction mapping...'
            % (time.time() - t))

        print('Writing output...')
        sys.stdout.flush()
        t = time.time()
        # Write out ROC curve results
        if actual_reaction_mapping is not None:
            with open(self._args.outpath +
                      '/roc/reaction_bayes.tsv', 'w') as f:
                write_roc_curve(f, rxn_bayes_pred.model1.reactions,
                                rxn_bayes_pred.model2.reactions,
                                rxn_bayes_pred,
                                actual_reaction_mapping)

        # Parse and output raw mapping
        if self._args.raw:
            rxn_bayes_pred.get_raw_map().to_csv(
                self._args.outpath + '/bayes_reactions.tsv', sep='\t',
                encoding='utf-8')

        # Output best mapping
        reaction_best = rxn_bayes_pred.get_best_map(
            self._args.threshold_reaction)
        reaction_best.to_csv(
            self._args.outpath + '/bayes_reactions_best.tsv', sep='\t',
            encoding='utf-8')
        print('It took %.2f seconds to write output...' % (time.time() - t))
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

        # read compartment map
        compartment_map = {}
        if self._args.compartment_map is not None:
            with open(self._args.compartment_map) as f:
                for row in f:
                    old, new = row.strip().split()
                    compartment_map[old] = new

        if self._args.compound_only:
            self._curate_compound(curator, model1, model2)
        else:
            self._curate_reaction(curator, model1, model2, compartment_map)

    def _curate_reaction(self, curator, model1, model2, compartment_map):
        print('Starting to do manual curation...\n')
        print('Type "stop" to stop\n')
        # iterate through reaction mapping
        for rmap in curator.reaction_map.iterrows():
            # skip already curated reactions
            if (curator.reaction_checked(rmap[0]) or
                    curator.reaction_ignored(rmap[0][0])):
                continue
            print('You have mapped %i reactions. '
                  'There are %i reactions unmapped.'
                  % (curator.num_curated_reactions,
                     curator.num_curated_reactions_left))
            print('\n')
            # check the compound mapping in current reaction
            compounds = curation.search_reaction(model1, [rmap[0][0]])
            compounds = [c.name for c in next(compounds)]
            dest_compounds = curation.search_reaction(model2, [rmap[0][1]])
            dest_compounds = [c.name for c in next(dest_compounds)]
            print('\n', '*' * 60, '\n')
            print('Below are the compound mapping involved:\n')
            for compound in compounds:
                for cmap in curator.compound_map.loc[compound].iterrows():
                    if (cmap is None or
                            cmap[0] not in dest_compounds or
                            curator.compound_checked((compound, cmap[0])) or
                            curator.compound_ignored(compound)):
                        continue
                    print('You have mapped %i compounds. '
                          'There are %i compounds unmapped.'
                          % (curator.num_curated_compounds,
                             curator.num_curated_compounds_left))
                    print('\n')
                    print('-' * 30)
                    print(cmap[1])
                    print('-' * 30)
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
                            curator.num_curated_compounds += 1
                            curator.num_curated_compounds_left -= 1
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
            print('You have mapped %i reactions. '
                  'There are %i reactions unmapped.'
                  % (curator.num_curated_reactions,
                     curator.num_curated_reactions_left))
            print('\n')
            print('Here is the reaction mapping:')
            print('-' * 30)
            print(rmap[1])
            print('-' * 30)
            print('\n')
            compounds = curation.search_reaction(model1, [rmap[0][0]])
            compounds = [c for c in next(compounds)]
            dest_compounds = curation.search_reaction(model2, [rmap[0][1]])
            dest_compounds = [c for c in next(dest_compounds)]
            print(('These two reactions have the following curated compound '
                   'pairs:\n'))
            count = 0
            compounds_curated = set()
            dest_compounds_curated = set()
            for compound in compounds:
                if compound.name in curator.curated_compound_map.index:
                    for cmap in curator.curated_compound_map.loc[
                            compound.name].iterrows():
                        compartment = compartment_map.get(
                            compound.compartment, compound.compartment
                        )
                        dest_compound = Compound(cmap[0], compartment)
                        if dest_compound in dest_compounds:
                            print(compound, dest_compound, cmap[1]['p'])
                            compounds_curated.add(compound)
                            dest_compounds_curated.add(dest_compound)
                            count += 1
            print('\n')
            print('%i compounds in %s' % (len(compounds), rmap[0][0]))
            unpaired = set(compounds) - compounds_curated
            if len(unpaired) > 0:
                print('Unpaired compounds: %s\n' % ', '.join(
                    [str(c) for c in unpaired]))
            print('%i compounds in %s' % (len(dest_compounds), rmap[0][1]))
            unpaired = set(dest_compounds) - dest_compounds_curated
            if len(unpaired) > 0:
                print('Unpaired compounds: %s\n' % ', '.join(
                    [str(c) for c in unpaired]))
            print('%i curated compound pairs\n' % count)
            ask = ''
            # waiting for legal curation input
            while ask not in ['y', 'n', 'stop']:
                ask = input(
                    ('True reaction match? (y/n/ignore/save/stop, '
                     'type ignore to ignore this reaction in future, '
                     'type save to save current progress, '
                     'type stop to save and exit): ')
                )
                if ask.lower() == 'n':
                    curator.add_mapping(rmap[0], 'r', False)
                    break
                if ask.lower() == 'y':
                    curator.add_mapping(rmap[0], 'r', True)
                    curator.num_curated_reactions += 1
                    curator.num_curated_reactions_left -= 1
                    break
                if ask.lower() == 'ignore':
                    curator.add_ignore(rmap[0][0], 'r')
                    break
                if ask.lower() == 'stop':
                    curator.save()
                    exit(0)
                if ask.lower() == 'save':
                    curator.save()
            print('\n', '=' * 60, '\n')
        curator.save()

    def _curate_compound(self, curator, model1, model2):
        for cmap in curator.compound_map.iterrows():
            if (curator.compound_checked(cmap[0]) or
                    curator.compound_ignored(cmap[0][0])):
                continue
            print('You have mapped %i compounds. '
                  'There are %i compounds unmapped.'
                  % (curator.num_curated_compounds,
                     curator.num_curated_compounds_left))
            print('-' * 30)
            print(cmap[1])
            print('-' * 30)
            print('\n')
            curation.search_compound(model1, [cmap[0][0]])
            curation.search_compound(model2, [cmap[0][1]])
            # waiting for legal curation input
            while True:
                ask = input(
                    ('True compound match? (y/n/ignore/save/stop, '
                     'type ignore to ignore this compound in future, '
                     'type save to save current progress, '
                     'type stop to save and exit): '))
                if ask.lower() == 'n':
                    curator.add_mapping(
                        cmap[0],
                        'c',
                        False
                    )
                    break
                if ask.lower() == 'y':
                    curator.add_mapping(
                        cmap[0],
                        'c',
                        True
                    )
                    curator.num_curated_compounds += 1
                    curator.num_curated_compounds_left -= 1
                    break
                if ask.lower() == 'ignore':
                    curator.add_ignore(cmap[0][0], 'c')
                    break
                if ask.lower() == 'stop':
                    curator.save()
                    exit(0)
                if ask.lower() == 'save':
                    curator.save()
            print('\n', '=' * 60, '\n')

    def _translate_id(self):
        """Translate ids"""
        # make output dir
        mkdir_p(self._args.outpath)

        # read mapping files
        cpd_mapping_id = tr_id.read_mapping(self._args.compound_map)
        rxn_mapping_id = tr_id.read_mapping(self._args.reaction_map)
        if self._args.compartment_map is not None:
            compartment_mapping_id = tr_id.read_mapping(
                self._args.compartment_map)
        else:
            compartment_mapping_id = None

        # make output model
        out_nm = tr_id.TranslatedModel(
            self._model,
            cpd_mapping_id,
            rxn_mapping_id,
            compartment_mapping_id
        )

        out_nm.write_model(dest=self._args.outpath, split_subsystem=False)


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
