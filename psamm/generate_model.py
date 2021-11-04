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
# Copyright 2019-2021  Christopher Powers <c-11060@uri.edu>

from __future__ import unicode_literals
import logging
import argparse
import yaml
import os
import sys
from collections import defaultdict
from collections import OrderedDict
import re
from six import iteritems
from Bio.KEGG import REST
from Bio.KEGG import Enzyme
from psamm.datasource.native import parse_reaction_equation_string
from psamm.datasource.native import NativeModel, ModelReader, ModelWriter
from psamm.expression import boolean
from psamm.datasource.entry import (DictReactionEntry as ReactionEntry)
from psamm.datasource.context import FileMark
from psamm.datasource.reaction import Reaction, Compound
from psamm.expression.affine import Expression
from psamm.formula import Formula
import time
logger = logging.getLogger(__name__)




class InputError(Exception):
    """Exception used to signal a general input error."""

def dict_representer(dumper, data):
    return dumper.represent_dict(data.items())

def dict_constructor(loader, node):
    return OrderedDict(loader.construct_pairs(node))

def parse_orthology(orthology, type):
    # Dictionary of reactions to genes
    asso_dict=defaultdict(lambda:[])
    # Populate the dictionary
    with open(orthology, "r") as infile:
        for line in infile:
            if line.startswith('#'):
                continue
            line=line.rstrip()
            listall=re.split("\t", line)
            if type == "R":
                keys = listall[14].split(',')
            elif type == "EC":
                keys = listall[10].split(',')
            elif type == "KO":
                keys = listall[11].split(',')
            if len(keys) == 0:
                continue
            for k in keys:
                asso_dict[k].append(listall[0])

    return(asso_dict)

'''
Function to sort the downloaded kegg object into a format
that is compatible with the psamm api for storage in
a reactions.yaml file.
'''
def model_reactions(reaction_entry_list):
        for reaction in reaction_entry_list:
            d = OrderedDict()
            d['id'] = encode_utf8(reaction.id)

            enzymes_list = []
            for i in reaction.enzymes:
                enzymes_list.append(i)

            pathways_list = []
            if reaction.pathways is not None:
                for i in reaction.pathways:
                    pathways_list.append(i[1])
            if len(pathways_list) == 0:
                pathways_list = None

            orth_list = []
            for i in reaction.orthology:
                    orth_list.append(i)
            if len(orth_list) == 0:
                orth_list = None

            if hasattr(reaction, 'name') and reaction.name is not None:
                d['name'] = encode_utf8(reaction.name)
            if hasattr(reaction, 'names') and reaction.names is not None:
                names_l = []
                for i in reaction.names:
                    names_l.append(i)
                d['names'] = encode_utf8(names_l)
            if hasattr(reaction, 'equation') and reaction.equation is not None:
                d['equation'] = encode_utf8(str(reaction.equation))
            if hasattr(reaction, 'enzymes') and reaction.enzymes is not None:
                d['enzymes'] = encode_utf8(enzymes_list)
            if hasattr(reaction, 'pathways') and reaction.pathways is not None:
                d['pathways'] = encode_utf8_list(pathways_list)
            if hasattr(reaction, 'comment') and reaction.comment is not None:
                d['comment'] = encode_utf8(str(reaction.comment))
            if hasattr(reaction, 'orthology') and reaction.orthology is not None:
                d['orthology'] = encode_utf8(orth_list)
            yield d

'''
This function handles specific edge cases in the kegg format of reaction
equations that are incompatible with the psamm format. Takes the downloaded
dictionary of reactions and returns the same dictionary with modified equations,
if necessary.
'''
def clean_reaction_equations(reaction_entry_list):
    for reaction in reaction_entry_list:
        equation = re.sub(r'\(.*?\)', lambda x: ''.join(x.group(0).split()), \
            str(reaction.equation))
        equation = equation.replace("(", "[")
        equation = equation.replace(")", "]")
        reaction.__dict__['values']['equation']=[equation]
    return(reaction_entry_list)

'''
This function creates a draft model based on the reactions:genes
file specified in the rxn variable. This function generates
reaction and compound information by utilizing the kegg
REST api to download the reaction information and uses psamm
functions to parse out the kegg data.
'''
def create_model_api(out, rxn_mapping):

    with open(os.path.join(out, "log.tsv"), "a+") as f:
        f.write("List of invalid Kegg IDs: \n")
    # Generate the yaml file for the reactions
    reaction_entry_list = []
    for entry in _download_kegg_entries(out, rxn_mapping, ReactionEntry):
        reaction_entry_list.append(entry)

    # Clean up the kegg reaction dict.
    reaction_entry_list = clean_reaction_equations(reaction_entry_list)

    # Sets up the yaml object for the reactions and writes
    # out the parsed reaction information to a reactions.yaml
    # file
    yaml.add_representer(OrderedDict, dict_representer)
    yaml.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
                         dict_constructor)
    yaml_args = {'default_flow_style': False,
                 'encoding': 'utf-8',
                 'allow_unicode': True}

    # Generate the yaml file for the compounds, filtering out the generic
    # compounds
    compound_entry_list=[]
    compound_set=set()
    for reaction in reaction_entry_list:
        eq = parse_reaction_equation_string(reaction.equation, 'c')
        for i in eq.compounds:
            compound_set.add(str(i[0].name))
    for entry in _download_kegg_entries(out, compound_set, CompoundEntry):
        compound_entry_list.append(entry)
    with open(os.path.join(out, 'compounds.yaml'), 'w+') as f:
        yaml.dump(list(model_compounds(compound_entry_list)), f, **yaml_args)

    # generate a yaml file for all generic compounds
    generic_compounds_list = check_generic_compounds(compound_entry_list)
    generic_entry_list=[]
    for entry in _download_kegg_entries(out, generic_compounds_list, CompoundEntry):
        generic_entry_list.append(entry)
    with open(os.path.join(out, 'compounds_generic.yaml'), 'w+') as f:
        yaml.dump(list(model_generic_compounds(generic_entry_list)), f, **yaml_args)
    with open(os.path.join(out, 'log.tsv'), "a+") as f:
        f.write("\nThere are {} generic compounds in the model\n".format(str(len(generic_compounds_list))))
        f.write("Generic compounds:\n")
        for i in generic_entry_list:
            f.write("{}".format(i.id))
            if i.name:
                f.write("|{}".format(i.name))
            else:
                f.write("|")
            if i.formula:
                f.write("|{}\n".format(i.formula))
            else:
                f.write("|\n")

    # generate a yaml for the reactions containing generic compounds and
    # and reactions without without generic compounds
    reaction_list_out = []
    reaction_list_generic = []
    with open(os.path.join(out, 'log.tsv'), 'a+') as f:
        f.write("\nThe reactions containing these generic compounds are: \n")
        for reaction in reaction_entry_list:
            if any(i in str(reaction.equation) for i in generic_compounds_list):
                f.write("{}".format(reaction.id))
                if reaction.name:
                    f.write("|{}".format(reaction.name))
                else:
                    f.write("|")
                if reaction.equation:
                    f.write("|{}\n".format(reaction.equation))
                else:
                    f.write("|\n")
                reaction_list_generic.append(reaction)
            else:
                reaction_list_out.append(reaction)
    with open(os.path.join(out, 'reactions.yaml'), 'w+') as f:
        yaml.dump(list(model_reactions(reaction_list_out)), f, **yaml_args)
    with open(os.path.join(out, 'reactions_generic.yaml'), 'w+') as f:
       	yaml.dump(list(model_reactions(reaction_list_generic)), f, **yaml_args)


    # Create a model.yaml file
    with open('{}/model.yaml'.format(out), mode='w') as myam:
        myam.write('default_flux_limit: 100\n')
        myam.write('default_compartment: {}\n'.format("c"))
        myam.write('extracellular: null\n')
        myam.write('biomass: null\n')
        myam.write('compounds:\n')
        myam.write('- include: ./compounds.yaml\n')
        myam.write('reactions:\n')
        myam.write('- include: ./reactions.yaml\n')
        myam.write('- include: ./gene-association.tsv\n')
        myam.write('  format: tsv\n')
        myam.write('model:\n')
        myam.write('- include: model_def.tsv\n')

    with open('{}/gene-association.tsv'.format(out), mode='w') as outfile:
        outfile.write('id\tgenes\n')
        for reaction in reaction_list_out:
            if len(rxn_mapping[reaction.id]) == 1:
                gene_asso = rxn_mapping[reaction.id]
            else:
                gene_asso = ['({})'.format(gene) for gene in rxn_mapping[reaction.id]]
            outfile.write('{}\t{}\n'.format(reaction.id, ' or '.join(gene_asso)))

    with open('{}/gene-association_generic.tsv'.format(out), mode='w') as outfile:
        outfile.write('id\tgenes\n')
        for reaction in reaction_list_generic:
            if len(rxn_mapping[reaction.id]) == 1:
                gene_asso = rxn_mapping[reaction.id]
            else:
                gene_asso = ['({})'.format(gene) for gene in rxn_mapping[reaction.id]]
            outfile.write('{}\t{}\n'.format(reaction.id, ' or '.join(gene_asso)))

    # Create a model_def file with all of the new reactions
    with open('{}/model_def.tsv'.format(out), mode='w') as f:
        for reaction in reaction_list_out:
            f.write('{}\n'.format(reaction.id))

    # Write out some final statistics for the model
    with open(os.path.join(out, 'log.tsv'), 'a+') as f:
        f.write("\nThere are {} reactions in the model".format(str(len(reaction_list_out))))
       	f.write("\nThere are {} compounds in the model\n".format(str(len(compound_entry_list))))




'''
Function to sort the downloaded	kegg object into a format
that is	compatible with the psamm api for storage in
a compounds.yaml file.
'''
def model_compounds(compound_entry_list):
    non_gen_compounds = []
    for compound in compound_entry_list:
        try:
            form = Formula.parse(str(compound.formula))
            if form.is_variable():
                continue
            elif compound.formula is None:
                continue
            elif 'R' in str(compound.formula):
                continue
            else:
                d = OrderedDict()
                d['id'] = encode_utf8(compound.id)
                non_gen_compounds.append(compound.id)
                if hasattr(compound, 'name') and compound.name is not None:
                    d['name'] = encode_utf8(compound.name)
                if hasattr(compound, 'names') and compound.names is not None:
                    names_l = []
                    for i in compound.names:
                        names_l.append(i)
                    d['names'] = encode_utf8(names_l)
                if hasattr(compound, 'formula') and compound.formula is not None:
                    d['formula'] = encode_utf8(str(compound.formula))
                if hasattr(compound, 'mol_weight') and compound.mol_weight is not None:
                    d['mol_weight'] = encode_utf8(compound.mol_weight)
                if hasattr(compound, 'comment') and compound.comment is not None:
                    d['comment'] = encode_utf8(str(compound.comment))
                if hasattr(compound, 'dblinks') and compound.dblinks is not None:
                    for key, value in compound.dblinks:
                        d['{}'.format(key)] = encode_utf8(value)
                yield d
        except:
            continue

def check_generic_compounds(compound_entry_list):
    generic_compounds_list = []
    for compound in compound_entry_list:
        try:
            form = Formula.parse(str(compound.formula))
            if form.is_variable():
                generic_compounds_list.append(compound.id)
            elif compound.formula is None:
                generic_compounds_list.append(compound.id)
            elif 'R' in str(compound.formula):
                generic_compounds_list.append(compound.id)
            else:
                continue
        except:
            generic_compounds_list.append(compound.id)
    return(generic_compounds_list)

def model_generic_compounds(compound_entry_list):
        #generic_compounds_list = []
        non_gen_compounds = []
        for compound in compound_entry_list:
            try:
                form = Formula.parse(str(compound.formula))
                d = OrderedDict()
                d['id'] = encode_utf8(compound.id)
                non_gen_compounds.append(compound.id)
                if hasattr(compound, 'name') and compound.name is not None:
                    d['name'] = encode_utf8(compound.name)
                if hasattr(compound, 'names') and compound.names is not None:
                    names_l = []
                    for i in compound.names:
                        names_l.append(i)
                    d['names'] = encode_utf8(names_l)
                if hasattr(compound, 'formula') and compound.formula is not None:
                    d['formula'] = encode_utf8(str(compound.formula))
                if hasattr(compound, 'mol_weight') and compound.mol_weight is not None:
                    d['mol_weight'] = encode_utf8(compound.mol_weight)
                if hasattr(compound, 'comment') and compound.comment is not None:
                    d['comment'] = encode_utf8(str(compound.comment))
                if hasattr(compound, 'dblinks') and compound.dblinks is not None:
                    for key, value in compound.dblinks:
                        d['{}'.format(key)] = encode_utf8(value)
                yield d
            except:
                print(compound)
                continue


'''
Two functions listed here for compatibility with
python 2 and python 3
'''
def encode_utf8(s):
    is_python2 = sys.version_info.major == 2
    if is_python2:
        if isinstance(s, unicode):
            return s.encode('utf-8')
    else:
        return s
def encode_utf8_list(s):
    is_python2 = sys.version_info.major == 2
    if is_python2:
        if isinstance(s, unicode):
            return s[1].encode('utf-8')
    else:
        return s

'''
Downloads the kegg entry associated with a reaction or
a compound and stores each line in an object that can
be parsed as a reaction or a compound, depending on the
input
'''
def _download_kegg_entries(out, rxn_mapping, entry_class, context=None):
    # go through the rxn mapping dict and pull all of the
    # reaction information out of the Kegg API
    with open(os.path.join(out,'log.tsv'), "a+") as f:
        for reactions in rxn_mapping:
            # Check for standard formatting of R and C. These can be treated
            # the same way
            if reactions[0]=="R" or reactions[0]=="C":
                # Try except loop to catch a failure to download.
                # On a failure, the failutr gets recorded in the log
                try:
                    request = REST.kegg_get(reactions)
                except:
                    f.write("".join(["  - ",reactions,"\n"]))
                    # sleep to give time to exit code
                    time.sleep(0.5)
                    continue
                entry_line = None
                section_id = None
                reaction = {}
                for lineno, line in enumerate(request):
                    line=line.rstrip()
                    if line=="///":
                        continue
                    if entry_line is None:
                        entry_line = lineno
                    # Look for the beginning of the section
                    m = re.match(r'([A-Z_]+)\s+(.*)', line.rstrip())
                    if m is not None:
                        section_id = m.group(1).lower()
                        reaction[section_id] = [m.group(2)]
                    elif section_id is not None:
                        reaction[section_id].append(line.strip())
                    else:
                        raise ParseError(
                            'Missing section identifier at line \
                            {}'.format(lineno))
                mark = FileMark(context, entry_line, 0)
                yield entry_class(reaction, filemark=mark)

                mark = FileMark(context, entry_line, 0)

"""
Functions converts gene associations to EC into gene associations for reaction
IDs. Returns a dictionary of Reaction IDs to genes.
"""
def parse_rxns_from_EC(rxn_mapping):
    rxn_dict=defaultdict(lambda:[])
    for reactions in rxn_mapping:
        try:
            request = REST.kegg_get(reactions)
        except:
            f.write("".join(["  - ",reactions,"\n"]))
            continue
        entry_line = None
        section_id = None
        reaction = {}
        for lineno, line in enumerate(request):
            line=line.rstrip()
            if line=="///":
                continue
            if entry_line is None:
                entry_line = lineno
            # Look for the beginning of the section
            m = re.match(r'([A-Z_]+)\s+(.*)', line.rstrip())
            if m is not None:
                section_id = m.group(1).lower()
                reaction[section_id] = [m.group(2)]
            elif section_id is not None:
                reaction[section_id].append(line.strip())
            else:
                raise ParseError(
                    'Missing section identifier at line \
                    {}'.format(lineno))
            rxn_ids=[]
        for i in reaction:
            if i == "all_reac":
                listall=re.split("\s",reaction[i][0])
                for r in listall:
                    if r[0]=="R":
                        rxn_dict[r]+=rxn_mapping[reactions]
    return(rxn_dict)


def main(args=None):
    """Entry point for the model generation script"""

    parser = argparse.ArgumentParser(
        description="Generate model based on Eggnog Annotations")
    parser.add_argument('--annotations', metavar='path',
        help='Path to the annotation file from Eggnog')
    parser.add_argument('--type', type=str,
        help='''Define whether to build the model on reaction ID, KO, or EC.\n
        options are: [R, KO, EC]''')
    parser.add_argument('--out', metavar='out',
        help='''Path to the output location for the model directory. This will\n
        be created for you.''')
    args = parser.parse_args(args)

    # Set up logging for the command line interface
    if 'PSAMM_DEBUG' in os.environ:
        level = getattr(logging, os.environ['PSAMM_DEBUG'].upper(), None)
        if level is not None:
            logging.basicConfig(level=level)
    else:
        logging.basicConfig(
            level=logging.INFO, format='%(levelname)s: %(message)s')


    # Check the validity of the input information
    if not args.annotations:
        raise InputError('Please specify a path to the eggnog annotations')
    if not args.type:
        raise InputError('Please specify one of R, KO, or EC as the type')
    if not args.out:
        raise InputError('Please specify an output directory')
    if os.path.exists(args.out):
        raise InputError('The output directory already exists! Exiting...')
        exit()
    else:
        os.mkdir(args.out)

    # Check the format of the eggnog annotation file.
    ## Add some code here to check the input file.

    # Launch the parse_orthology function to contruct a gene association file
    ortho_dict = parse_orthology(args.annotations, args.type)

    # convert EC to reactions
    if args.type == "EC":
        ortho_dict = parse_rxns_from_EC(ortho_dict)

    # Create the model using the kegg api
    create_model_api(args.out, ortho_dict)


class ReactionEntry(object):
    """Representation of entry in KEGG reaction file"""

    def __init__(self, values, filemark=None):
        self.values = dict(values)
        if 'entry' not in values:
            raise ParseError('Missing reaction identifier')
        self._id, _ = values['entry'][0].split(None, 1)
        self._filemark = filemark

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        try:
            return next(self.names)
        except StopIteration:
            return None

    @property
    def names(self):
        if 'name' in self.values:
            for line in self.values['name']:
                for name in line.rstrip(';').split(';'):
                    yield name.strip()

    @property
    def definition(self):
        if 'definition' not in self.values:
            return None
        return self.values['definition'][0]

    @property
    def equation(self):
        if 'equation' not in self.values:
            return None
        return self.values['equation'][0]

    @property
    def enzymes(self):
        if 'enzyme' in self.values:
            for line in self.values['enzyme']:
                for enzyme in line.split():
                    yield enzyme

    @property
    def pathways(self):
        if 'pathway' in self.values:
            for line in self.values['pathway']:
                pathway, name = line.split(None, 1)
                yield pathway, name

    @property
    def comment(self):
        if 'comment' not in self.values:
            return None
        return '\n'.join(self.values['comment'])

    @property
    def orthology(self):
        if 'orthology' in self.values:
            KO_name_dict = {}
            for line in self.values['orthology']:
                split = line.split(None, 1)
                yield split[0]

    def __getitem__(self, name):
        if name not in self.values:
            raise AttributeError('Attribute does not exist: {}'.format(name))
        return self._values[name]

    @property
    def filemark(self):
        return self._filemark

    def __repr__(self):
        return '<ReactionEntry "{}">'.format(self.id)

class ParseError(Exception):
    """Exception used to signal errors while parsing"""

class CompoundEntry(object):
    """Representation of entry in KEGG compound file"""

    def __init__(self, values, filemark=None):
        self.values = dict(values)
        if 'entry' not in values:
            raise ParseError('Missing compound identifier')
        self._id, _ = values['entry'][0].split(None, 1)
        self._filemark = filemark

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        try:
            return next(self.names)
        except StopIteration:
            return None

    @property
    def names(self):
        if 'name' in self.values:
            for line in self.values['name']:
                for name in line.rstrip(';').split(';'):
                    yield name.strip()

    @property
    def reactions(self):
        if 'reaction' in self.values:
            for line in self.values['reaction']:
                for rxnid in line.split():
                    yield rxnid

    @property
    def enzymes(self):
        if 'enzyme' in self.values:
            for line in self.values['enzyme']:
                for enzyme in line.split():
                    yield enzyme

    @property
    def formula(self):
        if 'formula' not in self.values:
            return None
        return self.values['formula'][0]

    @property
    def mol_weight(self):
        if 'mol_weight' not in self.values:
            return None
        return float(self.values['mol_weight'][0])

    @property
    def comment(self):
        if 'comment' not in self.values:
            return None
        return '\n'.join(self.values['comment'])

    def __getitem__(self, name):
        if name not in self.values:
            raise AttributeError('Attribute does not exist: {}'.format(name))
        return self._values[name]

    @property
    def filemark(self):
        return self._filemark

    def __repr__(self):
        return '<CompoundEntry "{}">'.format(self.id)
