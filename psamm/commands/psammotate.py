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
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>
# Copyright 2019-2020 Jing Wang <jingwang89@uri.edu>

from __future__ import print_function
from collections import defaultdict
import argparse
import re
import csv
import logging
from os import path, mkdir

from ..command import Command
from ..expression import boolean
from psamm.datasource.native import ModelWriter
from psamm.importer import write_yaml_model
logger = logging.getLogger(__name__)


class PsammotateCommand(Command):
    """Generate draft model from a template using provided gene mapping.

    This function will take a provided reciprocal best hit file
    """

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--rbh', metavar='file', type=argparse.FileType('r'),
            help=('A two column mapping of genes from the template '
                  'organism to the target organism.'))
        parser.add_argument(
            '--template', type=int,
            help=('The column of the RBH file where the template model '
                  'genes are listed, starts from 1.'))
        parser.add_argument(
            '--target', type=int,
            help=('The column of the RBH file where the target '
                  'model genes are listed, starts from 1.'))
        parser.add_argument(
            '--ignore-na', action='store_true',
            help=('Ignore the reactions that do not have gene association. '
                  '(default: retain these reactions in new reactions file)'))
        parser.add_argument(
            '--suffix', type=str, default=None,
            help='Suffix to append to end of reaction IDs and compartments.')
        g = parser.add_mutually_exclusive_group()
        g.add_argument(
            '--export-model', type=str, default=None,
            help='Path to directory for full model export. Cannot be specified'
                 ' with --output option.')
        g.add_argument(
            '--output', type=str, default='homolo_reactions',
            help=('The prefix of output YAML file, '
                  '(default: homolo_reactions). Cannot be specified with '
                  '--export-model option.'))
        super(PsammotateCommand, cls).init_parser(parser)

    def run(self):
        """Run psammotate command"""

        if self._args.export_model is None:
            if path.exists('{}.yaml'.format(self._args.output)):
                logger.warning('File {}.yaml already exists. '
                               'Please choose a different file name '
                               'through the --output '
                               'option.'.format(self._args.output))
                quit()
        elif self._args.export_model is not None:
            if path.exists('{}'.format(self._args.export_model)):
                logger.warning('Output directory {} already exists.'.format(
                    self._args.export_model))
                quit()
        if self._args.suffix:
            if any(i in self._args.suffix for i in [' ', '|', ':', ';', ',']):
                logger.error('Special character or space found '
                             'in provided suffix')
        tdict = app_reader(self._args.rbh, self._args.target - 1,
                           self._args.template - 1)
        trans_genes_dict = model_loader(self._model,
                                        self._args.ignore_na, tdict,
                                        self._args.suffix)
        homolo_reactions = list()
        for r in self._model.reactions:
            if r.id in self._model.model:
                if trans_genes_dict[r][2]:
                    r.properties['genes'] = remove_gap(
                        str(trans_genes_dict[r][1]))
                    if self._args.suffix:
                        r.properties['id'] += '_{}'.format(self._args.suffix)
                        for cpd, v in r.equation.compounds:
                            cpd._compartment += '_{}'.format(self._args.suffix)
                    homolo_reactions.append(r)
        if self._args.export_model is None:
            with open(self._args.output + '.yaml', 'w') as o:
                ModelWriter().write_reactions(o, homolo_reactions)
        elif self._args.export_model is not None:
            mkdir('{}'.format(self._args.export_model))
            reaction_list = [i.id for i in homolo_reactions]
            original_reactions = [i.id for i in self._model.reactions]
            for reaction in original_reactions:
                self._model.reactions.discard(reaction)
            for reaction in homolo_reactions:
                self._model.reactions.add_entry(reaction)
            compound_set = set()
            for reaction in self._model.reactions:
                for compound in reaction.equation.compounds:
                    compound_set.add(compound[0].name)
            compound_removal_set = set()
            for compound in self._model.compounds:
                if compound.id not in compound_set:
                    compound_removal_set.add(compound.id)
            for cpd in compound_removal_set:
                self._model.compounds.discard(cpd)
            exchange_del_set = set()
            for key in self._model.exchange:
                if key.name not in compound_set:
                    exchange_del_set.add(key)
            for ex in exchange_del_set:
                del self._model.exchange[ex]
            lim_del_set = set()
            for key in self._model.limits:
                if key not in reaction_list:
                    lim_del_set.add(key)
            for lim in lim_del_set:
                del self._model.exchange[lim]
            write_yaml_model(self._model, dest='{}'.format(
                self._args.export_model), split_subsystem=False)


def app_reader(app_file, query, template):
    '''This function will read in a gene mapping file and produce a
    gene to gene dictionary.

    The input is the app_file argument that the user enters, the query genome,
    which is the number of the genome that the user wishes to make a model of,
    and the number of the template genome which the user wishes to use to
    create the new model from.
    '''
    trans_dict = defaultdict(list)
    # Check what csv.reader in Python 3 takes (byte string or unicode string)
    for x, row in enumerate(csv.reader(app_file, delimiter='\t')):
        temp_l = []
        quer_l = []
        temp = row[template]
        quer = row[query]
        if temp != '-':
            if ',' in temp:
                temp_split = temp.split(',')
            else:
                temp_split = [temp]
            for i in temp_split:
                temp_l.append(i.strip())
            if ',' in quer:
                quer_split = quer.split(',')
            else:
                quer_split = [quer]
            for i in quer_split:
                quer_l.append(i.strip())
        for i in temp_l:
            for j in quer_l:
                trans_dict[i].append(j)
    return trans_dict


def model_loader(nm, ignore_na, translation_dict, suffix=None):
    '''This function will translate genes in one model
    based on a gene to gene mapping dictionary.

    The function takes a native model object, a true or false
    ignore_na argument, and a gene mapping dictionary translation_dict
    that is made by the app_reader function. The ignore_na argument
    will determine if reactions with no genes will be kept in the final
    translated model or not. This function will then translate the gene
    associations based on the gene mapping dictionary and produce a
    new set of gene associations based on the mapped associations. The
    gene expressions are evaluated based on the logic to see if they should
    be included in the new model.
    '''
    new_translation_dict = {}
    translated_genes = {}
    target_genes_l = {}
    target_model_reactions = []
    translation_dict.pop('-', None)
    for key, value in translation_dict.items():
        if len(value) > 1:
            value_s = '({})'.format(' or '.join(value))
        else:
            value_s, = value
        new_translation_dict[key] = value_s
        for i in value:
            target_genes_l[i] = True
    mm = nm.create_metabolic_model()
    model_rxns = [i for i in mm.reactions]
    print('ReactionID\tOriginal_Genes\tTranslated_Genes\tIn_final_model')
    for entry in nm.reactions:
        if entry.id in model_rxns:
            if entry.genes is None:
                target_model_reactions.append(entry.id)
                translated_genes[entry] = [None, None,
                                           not ignore_na]
                id = entry.id
                if suffix is not None:
                    id = id + '_{}'.format(suffix)
                print(u'{}\t{}\t{}\t{}'.format(
                    id, entry.genes, 'None',
                    not ignore_na))
            elif entry.genes is not None:
                genes = re.sub(r'\?', '', entry.genes)
                e = boolean.Expression(genes)
                e_1 = e.substitute(
                    lambda v: new_translation_dict.get(
                        v.symbol, '-') != '-')
                genes_1 = entry.genes
                gene_list = get_gene_list(genes)
                for g in gene_list:
                    if g in new_translation_dict:
                        genes = re.sub(r'\b' + g + r'\b',
                                       new_translation_dict[g], genes)
                    else:
                        genes = re.sub(r'\b' + g + r'\b',
                                       '-', genes)
                id = entry.id
                if suffix is not None:
                    id = id + '_{}'.format(suffix)
                print(u'{}\t{}\t{}\t{}'.format(
                    id, entry.genes, genes, e_1.value))
                translated_genes[entry] = [genes_1, genes, e_1.value]
    return translated_genes


def get_gene_list(genes):
    '''This function is used to get a list of genes from a gene expression
    '''
    gene_list = re.sub(r'[,()\']*', '', genes)
    gene_list = re.sub(r'\band\b', '', gene_list)
    gene_list = re.sub(r'\bor\b', '', gene_list)
    gene_list = re.split(r'\s+', gene_list)  # Transfer string to list
    gene_list = frozenset(gene_list)
    return gene_list


def remove_gap(string):
    '''This function is used to simplify translated gene expression strings.
    '''
    while True:
        prev = string
        string = re.sub(r'or - or', 'or', string)
        string = re.sub(r' or -', '', string)
        string = re.sub(r'\(-\)', '-', string)
        string = re.sub(r'- or ', '', string)
        string = re.sub(r'\(\)', '', string)
        string = re.sub(r'\'\'', '', string)
        if prev == string:
            return string
