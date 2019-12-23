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
# Copyright 2014-2017  Keith Dufault-Thompson <keitht547@uri.edu>
# Copyright 2019    Jing Wang <jingwang89@uri.edu>

from __future__ import print_function
from collections import defaultdict, OrderedDict
import argparse
import re
import csv
import logging
import os
import yaml

from ..command import Command
from ..expression import boolean
from psamm.datasource.native import ModelWriter
logger = logging.getLogger(__name__)


class PsammotateCommand(Command):
    """Generate draft model from a template using provided gene mapping.

    This function will take a provided reciprocal best hit file
    """

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--rbh', metavar='file', type=argparse.FileType('r'),
            help=('The homologous gene table file, generally '
                  'from pan-genome analysis'))
        parser.add_argument(
            '--template', type=int,
            help=('The column of the RBH file where the template model '
                  'genes are listed, starts from 1'))
        parser.add_argument(
            '--target', type=int,
            help=('The column of the RBH file where the target '
                  'model genes are listed, starts from 1'))
        parser.add_argument(
            '--prefix', type=str, default='homolo_reactions',
            help=('The prefix of output YAML file, '
                  '(default: homolo_reactions)'))
        parser.add_argument(
            '--ignore-na', action='store_true',
            help=('Ignore the reactions that do not have gene association. '
                  '(default: retain these reactions in new model)'))
        super(PsammotateCommand, cls).init_parser(parser)

    def run(self):
        """Run psammotate command"""
        # print(self._args.rbh)
        for i in app_reader(self._args.rbh, self._args.target - 1,
                            self._args.template - 1):
            tdict = i
        for i in model_loader(self, tdict):
            trans_genes_dict = i
        # print_yaml(self, trans_genes_dict)
        homolo_reactions = list()
        for r in self._model.reactions:
            if trans_genes_dict[r][2]:
                r.properties['genes'] = trans_genes_dict[r][1]
                homolo_reactions.append(r)
        with open(self._args.prefix + '.yaml', 'w') as o:
            ModelWriter().write_reactions(o, homolo_reactions)


def app_reader(app_file, query, template):
    '''This function will read in the app file and produce the mapping
    dictionary.

    The input is the app_file argument that the user enters, the query genome,
    which is the number of the genome that the user wishes to make a model of,
    and the number of the template genome which the user wishes to use to
    create the new model from.
    '''

    new_list = {}
    trans_dict = defaultdict(list)
    transl_dict = {}
    # Check what csv.reader in Python 3 takes (byte string or unicode string)
    for x, row in enumerate(csv.reader(app_file, delimiter='\t')):
        temp_l = []
        quer_l = []
        if x == 0:
            new_list['template'] = (row[template])
            new_list['query'] = (row[query])
        else:

            temp = row[template]
            quer = row[query]
            if temp != '-':
                temp_split = []
                # print(temp, quer)
                if ',' in temp:
                    temp_split = temp.split(',')
                else:
                    temp_split.append(temp)
                for i in temp_split:
                    if ':' in i:
                        i = i.split(':')[1]
                    else:
                        i = i
                    temp_l.append(i)
                if ',' in quer:
                    quer_split = quer.split(',')
                else:
                    quer_split = [quer]
                for i in quer_split:
                    if ':' in i:
                        i = i.split(':')[1]
                    else:
                        i = i
                    quer_l.append(i)
        for i in temp_l:
            i = i.strip()
            for j in quer_l:
                j = j.strip()
                trans_dict[i].append(j)
        for key, value in trans_dict.items():
            value = str(value)
            value = value.strip('[')
            value = value.strip(']')
            value = value.strip()
            value = value.strip(',')
            value = value.replace(',', ' or')
            transl_dict[key] = value
    yield(trans_dict)


def model_loader(self, translation_dict):
    new_translation_dict = {}
    translated_genes = {}
    model = self._model
    target_genes_l = {}
    target_model_reactions = []
    translation_dict.pop('-', None)
    for key, value in translation_dict.items():
        value_s = None
        for i in value:
            if value_s is None:
                value_s = i
            elif value_s is not None:
                value_s = value_s + ' or ' + i
        if len(value) >= 2:
            value_s = '(' + value_s + ')'
        new_translation_dict[key] = value_s

    for value in translation_dict.values():
        for i in value:
            target_genes_l[i] = True
    target_genes_l.pop('-', None)
    target_genes_l.pop('\'-\'', None)
    target_genes_l['\'-\''] = False
    target_genes_l['-'] = False
    target_genes_l[None] = False
    target_genes_l['Sink'] = False
    target_genes_l['Gap'] = False
    target_genes_l['Synthesis'] = False
    target_genes_l['Diffusion'] = False
    target_genes_l['Biomass'] = False
    target_genes_l['()'] = False
    target_genes_l['None'] = False
    model_rxns = []
    mm = model.create_metabolic_model()
    for rxd in mm.reactions:
        model_rxns.append(rxd)
    for entry in model.reactions:
        if entry.id in model_rxns:
            # print(entry.id)
            # print('wp3_genes', entry.genes)
            if entry.genes is None:
                target_model_reactions.append(entry.id)
                translated_genes[entry] = [None, None,
                                           not self._args.ignore_na]
                print('{}\t{}\t{}\t{}\t{}'.format(
                    entry.id, entry.genes, 'None',
                    not self._args.ignore_na, not self._args.ignore_na))

            if entry.genes is not None:
                genes = re.sub(r'\?', '', entry.genes)
                e = boolean.Expression(genes)
                e_1 = e.substitute(
                    lambda v: new_translation_dict.get(
                        v.symbol, '-') != '-')
                # print('%s\t%s' % (entry.id, entry.genes))
                genes_1 = entry.genes
                gene_list = get_gene_list(genes)
                # print(gene_list)
                for g in gene_list:
                    if g in new_translation_dict:
                        genes = re.sub(r'\b' + g + r'\b',
                                       new_translation_dict[g], genes)
                    else:
                        genes = re.sub(r'\b' + g + r'\b',
                                       '-', genes)
                # genes = remove_gap(genes)
                print('{}\t{}\t{}\t{}\t{}'.format(
                    entry.id, entry.genes, genes, e_1, e_1.value))
                translated_genes[entry] = [genes_1, genes, e_1.value]
    yield(translated_genes)


def get_gene_list(genes):
    gene_list = re.sub(r'[,()\']*', '', genes)
    gene_list = re.sub(r'\band\b', '', gene_list)
    gene_list = re.sub(r'\bor\b', '', gene_list)
    gene_list = re.split(r'\s+', gene_list)  # Transfer string to list
    gene_list = frozenset(gene_list)
    return gene_list


def remove_gap(string):
    while True:
        prev = string
        string = re.sub(r'or - or', 'or', string)
        string = re.sub(r' or -', '', string)
        string = re.sub(r'- or ', '', string)
        string = re.sub(r'\(\)', '', string)
        string = re.sub(r'\'\'', '', string)
        if prev == string:
            return string


def encode_utf8(s):
    if isinstance(s, str):
        return s.encode('utf-8')
    return s


def print_yaml(self, trans_genes):
    yaml.add_representer(OrderedDict, dict_representer)

    def unicode_representer(dumper, uni):
        node = yaml.ScalarNode(tag=u'tag:yaml.org,2002:str', value=uni)
        return node
    yaml.add_representer(str, unicode_representer)
    yaml.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
                         dict_constructor)
    dest = '.'
    yaml_args = {'default_flow_style': False,
                 'encoding': 'utf-8',
                 'allow_unicode': True}
    with open(os.path.join(dest, 'homolo_reactions.yaml'), 'w+') as f:
        yaml.dump(list(model_export(self, trans_genes)), f, **yaml_args)


def dict_representer(dumper, data):
    return dumper.represent_dict(data.items())


def dict_constructor(loader, node):
    return OrderedDict(loader.construct_pairs(node))


def model_export(self, trans_genes):
    for entry, gene_list in trans_genes.items():
        if gene_list[2] is not False:
            # print(entry.id)

            d = OrderedDict()
            d['id'] = str(entry.id)

            if hasattr(entry, 'name') and entry.name is not None:
                d['name'] = encode_utf8(entry.name)
            if hasattr(entry, 'genes') and entry.genes is not None:
                d['genes'] = encode_utf8(gene_list[1])
            if hasattr(entry, 'equation') and entry.equation is not None:
                d['equation'] = encode_utf8(str(entry.equation))
            if hasattr(entry, 'subsystem') and entry.subsystem is not None:
                d['subsystem'] = encode_utf8(entry.subsystem)
            if hasattr(entry, 'ec') and entry.ec is not None:
                d['ec'] = encode_utf8(entry.ec)
            if hasattr(entry, 'subset') and entry.subset is not None:
                d['subset'] = encode_utf8(entry.subset)
            yield(d)
