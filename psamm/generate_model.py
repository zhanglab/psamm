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
# Copyright 2019-2022  Christopher Powers <c-11060@uri.edu>
# Copyright 2021-2022  Jason Vailionis <jason_vailionis@uri.edu>

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
from psamm.datasource.native import parse_reaction_equation_string
from psamm.datasource.native import ModelReader
from psamm.datasource.entry import (DictReactionEntry as ReactionEntry)
from psamm.datasource.context import FileMark
from psamm.datasource.reaction import Reaction, Compound, Direction
from psamm.formula import Formula, ParseError
from pkg_resources import resource_filename
import pkg_resources
from psamm.command import _trim
from six import add_metaclass, iteritems, text_type
import abc
from urllib.error import HTTPError
logger = logging.getLogger(__name__)
if sys.version_info.minor > 5:
    try:
        from Bio.KEGG import REST
    except ImportError:
        logger.warning("WARNING: Biopython package not found! "
                       "Some functions will be unusable")
    try:
        from libchebipy._chebi_entity import ChebiEntity
    except ImportError:
        logger.warning("WARNING: The Chebi API package not found! "
                       "Some functions will be unusable")


class ParseError2(Exception):
    """Exception used to signal errors while parsing"""


class VersionError(Exception):
    """Exception used to signal incorrect python version"""


class CommandError(Exception):
    """Error from running a command.

    This should be raised from a ``Command.run()`` if any arguments are
    misspecified. When the command is run and the ``CommandError`` is raised,
    the caller will exit with an error code and print appropriate usage
    information.
    """


class InputError(Exception):
    """Exception used to signal a general input error. Exiting..."""


def dict_representer(dumper, data):
    return dumper.represent_dict(data.items())


def dict_constructor(loader, node):
    return OrderedDict(loader.construct_pairs(node))


@add_metaclass(abc.ABCMeta)
class Command(object):
    """Represents a command in the interface, operating on a model.

    The constructor will be given the NativeModel and the command
    linenamespace. The subclass must implement :meth:`run` to
    handle commandexecution. The doc string will be used as
    documentation for the command in the command line interface.

    In addition, :meth:`init_parser` can be implemented as a
    classmethod which will allow the command to initialize an
    instance of :class:`argparse.ArgumentParser` as desired.
    The resulting argumentnamespace will be passed to the
    constructor.
    """

    def __init__(self, args):
        self._args = args

    @classmethod
    def init_parser(cls, parser):
        """Initialize command line parser (:class:`argparse.ArgumentParser`)"""

    @abc.abstractmethod
    def run(self):
        """Execute command"""

    def argument_error(self, msg):
        """Raise error indicating error parsing an argument."""
        raise CommandError(msg)

    def fail(self, msg, exc=None):
        """Exit command as a result of a failure."""
        logger.error(msg)
        if exc is not None:
            logger.debug('Command failure caused by exception!', exc_info=exc)
        sys.exit(1)


def parse_orthology(orthology, type, col):
    '''
    function to parse orthology tables. The default format for these
    tables is the default eggnog output, but custom colummn numers can be
    passed through the --col argument
    '''
    # Dictionary of reactions to genes
    asso_dict = defaultdict(lambda: [])
    # Populate the dictionary
    with open(orthology, "r") as infile:
        for line in infile:
            if line.startswith('#'):
                continue
            line = line.rstrip()
            listall = re.split("\t", line)
            # Add handling for a provided column number.
            if col:
                if len(listall[col-1]) > 0:
                    keys = listall[col-1].split(',')
                else:
                    keys = []
            else:
                if type == "R":
                    keys = listall[14].split(',')
                elif type == "EC":
                    keys = listall[10].split(',')
                elif type == "KO":
                    keys = listall[11].split(',')
                elif type == "tcdb":
                    if listall[17] == '-':
                        continue
                    keys = listall[17].split(',')
            if len(keys) == 0:
                continue
            for k in keys:
                asso_dict[k].append(listall[0])
    return(asso_dict)


def model_reactions(reaction_entry_list):
    '''
    Function to sort the downloaded kegg object into a format
    that is compatible with the psamm api for storage in
    a reactions.yaml file.
    '''
    for reaction in reaction_entry_list:

        d = OrderedDict()
        d['id'] = reaction.id

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
            d['name'] = reaction.name
        if hasattr(reaction, 'names') and reaction.names is not None:
            names_l = []
            for i in reaction.names:
                names_l.append(i)
            d['names'] = names_l
        if hasattr(reaction, 'equation') and reaction.equation is not None:
            d['equation'] = str(reaction.equation)
        if hasattr(reaction, 'definition') and reaction.definition is not None:
            d['KEGG_definition'] = reaction.definition
        if hasattr(reaction, 'enzymes') and reaction.enzymes is not None:
            d['enzymes'] = enzymes_list
        if hasattr(reaction, 'pathways') and reaction.pathways is not None:
            d['pathways'] = pathways_list
        if hasattr(reaction, 'comment') and reaction.comment is not None:
            d['comment'] = str(reaction.comment)
        if hasattr(reaction, 'tcdb_family') and \
                reaction.tcdb_family is not None:
            d['tcdb_family'] = str(reaction.tcdb_family)
        if hasattr(reaction, 'substrates') and \
                reaction.substrates is not None:
            d['substrates'] = reaction.substrates
        if hasattr(reaction, 'genes') and \
                reaction.genes is not None:
            d['genes'] = str(reaction.genes)
        if hasattr(reaction, 'orthology') and \
                reaction.orthology is not None:
            d['orthology'] = orth_list
        yield d


def clean_reaction_equations(reaction_entry_list):
    '''
    This function handles specific edge cases in the kegg
    format of reaction equations that are incompatible with
    the psamm format. Takes the downloaded dictionary of reactions
    and returns the same dictionary with modified equations,
    if necessary.
    '''
    generic = []
    for reaction in reaction_entry_list:
        equation = re.sub(r'\(.*?\)', lambda x: ''.join(x.group(0).split()),
                          str(reaction.equation))
        equation_out = []
        comp = re.split(" ", equation)
        comp = comp + [""]

        for i in range(0, len(comp)-1):
            # Handles unusual stoichiometry at the beginning of compounds
            if "(" in comp[i] and "C" in comp[i+1]:
                if i > 0:
                    if "=" in comp[i-1] or "+" in comp[i-1]:
                        equation_out.append(comp[i])
                        continue
                else:
                    equation_out.append(comp[i])
                    continue
                generic.append(reaction.id)
            # special handling to retain stoichiometry of (n+1) and (n-1), etc.
            # at the end of the compound
            if "(" in comp[i]:  # and (not "(n+" in i or not "(n-" in i):
                comp[i] = comp[i].replace("(", "[")
                comp[i] = comp[i].replace(")", "]")
                generic.append(reaction.id)
            equation_out.append(comp[i])
        equation = ' '.join(equation_out)
        reaction.__dict__['values']['equation'] = [equation]

    return(reaction_entry_list, generic)


def create_model_api(out, rxn_mapping, verbose, use_rhea, default_compartment):
    '''
    This function creates a draft model based on the reactions:genes
    file specified in the rxn variable. This function generates
    reaction and compound information by utilizing the kegg
    REST api to download the reaction information and uses psamm
    functions to parse out the kegg data.
    '''
    with open(os.path.join(out, "log.tsv"), "a+") as f:
        f.write("List of invalid Kegg IDs: \n")

    if verbose:
        logger.info("There are {} reactions to download"
                    ".".format(len(rxn_mapping)))
    # Generate the yaml file for the reactions
    reaction_entry_list = []
    count = 0
    for entry in _download_kegg_entries(out, rxn_mapping, None, ReactionEntry):
        reaction_entry_list.append(entry)
        count += 1
        if verbose and count % 25 == 0:
            logger.info("{}/{} reactions downloaded..."
                        .format(count, len(rxn_mapping)))

    # Clean up the kegg reaction dict.
    reaction_entry_list, gen = clean_reaction_equations(reaction_entry_list)

    # Sets up the yaml object for the reactions and writes
    # out the parsed reaction information to a reactions.yaml
    # file
    yaml.add_representer(OrderedDict, dict_representer)
    yaml.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
                         dict_constructor)
    yaml_args = {'default_flow_style': False,
                 'encoding': 'utf-8',
                 'allow_unicode': True}

    # Load database once if --rhea is given
    if use_rhea:
        rhea_db = RheaDb(resource_filename('psamm',
                         'external-data/chebi_pH7_3_mapping.tsv'))
    else:
        rhea_db = None

    compound_entry_list = []
    compound_set = set()
    for reaction in reaction_entry_list:
        eq = parse_reaction_equation_string(reaction.equation, 'c')
        for i in eq.compounds:
            compound_set.add(str(i[0].name))
    if verbose:
        logger.info("There are {} compounds to download."
                    .format(len(compound_set)))
    count = 0
    logger.info("Downloading kegg entries and performing charge correction")
    for entry in _download_kegg_entries(out, compound_set,
                                        rhea_db, CompoundEntry):
        compound_entry_list.append(entry)
        count += 1
        if verbose and count % 25 == 0:
            logger.info("{}/{} compounds downloaded...".
                        format(count, len(compound_set)))
    with open(os.path.join(out, 'compounds.yaml'), 'w+') as f:
        yaml.dump(list(model_compounds(compound_entry_list)), f, **yaml_args)

    # generate a yaml file for all generic compounds
    generic_compounds_list = check_generic_compounds(compound_entry_list)
    generic_entry_list = []
    logger.info("Downloading kegg entries for generic compounds and "
                "performing charge correction")
    for entry in _download_kegg_entries(out, generic_compounds_list, rhea_db,
                                        CompoundEntry):
        generic_entry_list.append(entry)
    with open(os.path.join(out, 'compounds_generic.yaml'), 'w+') as f:
        yaml.dump(list(model_generic_compounds(generic_entry_list)),
                  f, **yaml_args)
    with open(os.path.join(out, 'log.tsv'), "a+") as f:
        f.write("\nThere are {} generic compounds in the model\n".format(str(
            len(generic_compounds_list))))
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
    compound_list_out = set()

    with open(os.path.join(out, 'log.tsv'), 'a+') as f:
        f.write("\nThe reactions containing these generic compounds are: \n")
        for reaction in reaction_entry_list:
            if any(i in str(reaction.equation)
                   for i in generic_compounds_list) or reaction.id in gen:
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
                for c in str(reaction.equation):
                    compound_list_out.add(str(c))

    with open(os.path.join(out, 'reactions.yaml'), 'w+') as f:
        yaml.dump(list(model_reactions(reaction_list_out)), f, **yaml_args)
    with open(os.path.join(out, 'reactions_generic.yaml'), 'w+') as f:
        yaml.dump(list(model_reactions(reaction_list_generic)), f, **yaml_args)

    # Create a model.yaml file
    with open('{}/model.yaml'.format(out), mode='w') as myam:
        myam.write('default_flux_limit: 100\n')
        myam.write('default_compartment: {}\n'.format(default_compartment))
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
                gene_asso = ['({})'.format(gene)
                             for gene in rxn_mapping[reaction.id]]
            outfile.write('{}\t{}\n'.format(reaction.id,
                          ' or '.join(gene_asso)))

    with open('{}/gene-association_generic.tsv'.format(out), mode='w') \
            as outfile:
        outfile.write('id\tgenes\n')
        for reaction in reaction_list_generic:
            if len(rxn_mapping[reaction.id]) == 1:
                gene_asso = rxn_mapping[reaction.id]
            else:
                gene_asso = ['({})'.format(gene)
                             for gene in rxn_mapping[reaction.id]]
            outfile.write('{}\t{}\n'.format(reaction.id,
                          ' or '.join(gene_asso)))

    # Create a model_def file with all of the new reactions
    with open('{}/model_def.tsv'.format(out), mode='w') as f:
        for reaction in reaction_list_out:
            f.write('{}\n'.format(reaction.id))

    # Write out some final statistics for the model
    with open(os.path.join(out, 'log.tsv'), 'a+') as f:
        f.write("\nThere are {} reactions in the model".format(str(
                len(reaction_list_out))))
        f.write("\nThere are {} compounds in the model\n".format(str(
                len(compound_list_out))))
    if verbose:
        logger.info("\nThere are {} reactions in the model".format(str(
                    len(reaction_list_out))))
        logger.info("\nThere are {} compounds in the model\n".format(str(
                    len(compound_entry_list))))


def model_compounds(compound_entry_list):
    '''
    Function to sort the downloaded	kegg object into a format
    that is	compatible with the psamm api for storage in
    a compounds.yaml file.
    '''
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
                d['id'] = compound.id
                non_gen_compounds.append(compound.id)
                if hasattr(compound, 'name') and compound.name is not None:
                    d['name'] = compound.name
                if hasattr(compound, 'names') and \
                        compound.names is not None:
                    names_l = []
                    for i in compound.names:
                        names_l.append(i)
                    d['names'] = names_l
                if hasattr(compound, 'formula') and \
                        compound.formula is not None:
                    d['formula'] = str(compound.formula)
                if hasattr(compound, 'mol_weight') and \
                        compound.mol_weight is not None:
                    d['mol_weight'] = compound.mol_weight
                if hasattr(compound, 'comment') and \
                        compound.comment is not None:
                    d['comment'] = str(compound.comment)
                if hasattr(compound, 'dblinks') and \
                        compound.dblinks is not None:
                    for key, value in compound.dblinks:
                        if key != 'ChEBI':
                            d['{}'.format(key)] = value
                if hasattr(compound, 'chebi') \
                        and compound.chebi is not None:
                    d['ChEBI'] = compound.chebi
                if hasattr(compound, 'chebi_all') \
                        and compound.chebi_all is not None:
                    d['ChEBI_all'] = compound.chebi_all
                if hasattr(compound, 'charge') \
                        and compound.charge is not None:
                    d['charge'] = compound.charge
                yield d
        except ParseError:
            logger.warning("import of {} failed"
                           " and will not be imported into compounds.yaml "
                           "or compounds_generic.yaml".format(compound.id))
            continue


def check_generic_compounds(compound_entry_list):
    '''
    Function for checking if the compound formulation is
    compatible with psamm. generalized rules for this are
    that compounds must have a formula, the formula cannot
    be variable (e.g. presence of X), and R groups are
    generally discouraged.
    '''
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
        except ParseError:
            generic_compounds_list.append(compound.id)
    return(generic_compounds_list)


def model_generic_compounds(compound_entry_list):
    '''
    Function to sort the downloaded	kegg object into a format
    that is	compatible with the psamm api for storage in
    a generic_compounds.yaml file. This function contains
    special error handling for improperly formatted compounds
    '''
    non_gen_compounds = []
    for compound in compound_entry_list:
        try:
            d = OrderedDict()
            d['id'] = compound.id
            non_gen_compounds.append(compound.id)
            if hasattr(compound, 'name') and compound.name is not None:
                d['name'] = compound.name
            if hasattr(compound, 'names') and compound.names is not None:
                names_l = []
                for i in compound.names:
                    names_l.append(i)
                d['names'] = names_l
            if hasattr(compound, 'formula') and compound.formula is not None:
                d['formula'] = str(compound.formula)
            if hasattr(compound, 'mol_weight') and \
                    compound.mol_weight is not None:
                d['mol_weight'] = compound.mol_weight
            if hasattr(compound, 'comment') and compound.comment is not None:
                d['comment'] = str(compound.comment)
            if hasattr(compound, 'dblinks') and compound.dblinks is not None:
                for key, value in compound.dblinks:
                    if key != 'ChEBI':
                        d['{}'.format(key)] = value
            if hasattr(compound, 'chebi') and compound.chebi is not None:
                d['ChEBI'] = compound.chebi
            if hasattr(compound, 'chebi_all') and \
                    compound.chebi_all is not None:
                d['ChEBI_all'] = compound.chebi_all
            if hasattr(compound, 'charge') and compound.charge is not None:
                d['charge'] = compound.charge
            yield d
        except ParseError:
            logger.warning("{} is improperly formatted "
                           " and will not be imported into "
                           "compounds_generic.yaml".format(compound.id))
            continue


def _download_kegg_entries(out, rxn_mapping, rhea, entry_class, context=None):
    '''
    Downloads the kegg entry associated with a reaction or
    a compound and stores each line in an object that can
    be parsed as a reaction or a compound, depending on the
    input
    '''
    # go through the rxn mapping dict and pull all of the
    # reaction information out of the Kegg API
    with open(os.path.join(out, 'log.tsv'), "a+") as f:
        for reactions in rxn_mapping:
            # Check for standard formatting of R and C. These can be treated
            # the same way
            if reactions[0] == "R" or reactions[0] == "C":
                # Try except loop to catch a failure to download.
                # On a failure, the failure gets recorded in the log
                try:
                    request = REST.kegg_get(reactions)
                except HTTPError:
                    f.write("".join(["  - ", reactions, "\n"]))
                    continue
                entry_line = None
                section_id = None
                reaction = {}
                for lineno, line in enumerate(request):
                    line = line.rstrip()
                    if line == "///":
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
                        raise ParseError2(
                            'Missing section identifier at line \
                            {}'.format(lineno))
                mark = FileMark(context, entry_line, 0)
                yield entry_class(reaction, rhea, filemark=mark)


def parse_rxns_from_EC(rxn_mapping, out, verbose):
    """
    Functions converts gene associations to EC into gene
    associations for reaction IDs. Returns a dictionary
    of Reaction IDs to genes.
    """
    rxn_dict = defaultdict(lambda: [])
    if verbose:
        logger.info("Downloading reactions associated with EC...")
        logger.info("There are {} ECs download".format(len(rxn_mapping)))
        count = 0
    with open(os.path.join(out, "log.tsv"), "a+") as f:
        for reactions in rxn_mapping:
            if verbose:
                if count % 25 == 0:
                    logger.info("{}/{} have been downloaded".format(count,
                                len(rxn_mapping)))
                count += 1
            try:
                request = REST.kegg_get(reactions)
            except HTTPError:
                f.write("".join(["  - ", reactions, "\n"]))
                continue
            entry_line = None
            section_id = None
            reaction = {}
            for lineno, line in enumerate(request):
                line = line.rstrip()
                if line == "///":
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
                    raise ParseError2('Missing section identifier at line '
                                      '{}'.format(lineno))
            if "all_reac" in reaction:
                listall = re.split(" ", reaction["all_reac"][0])
                for r in listall:
                    if r[0] == "R":
                        rxn_dict[r] += rxn_mapping[reactions]

    return(rxn_dict)


def parse_rxns_from_KO(rxn_mapping, out, verbose):
    """
    Functions converts gene associations to EC into gene
    associations for reaction IDs. Returns a dictionary
    of Reaction IDs to genes.
    """
    rxn_dict = defaultdict(lambda: [])
    if verbose:
        logger.info("Downloading reactions associated with KO...")
        logger.info("There are {} KOs download".format(len(rxn_mapping)))
        count = 0
    with open(os.path.join(out, "log.tsv"), "a+") as f:
        for reactions in rxn_mapping:
            if verbose:
                if count % 25 == 0:
                    logger.info("{}/{} have been downloaded".format(count,
                                len(rxn_mapping)))
                count += 1
            try:
                request = REST.kegg_get(reactions)
            except HTTPError:
                f.write("".join(["  - ", reactions, "\n"]))
                continue
            entry_line = None
            section_id = None
            reaction = {}
            for lineno, line in enumerate(request):
                line = line.rstrip()
                if line == "///":
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
                    raise ParseError2(
                        'Missing section identifier at line \
                        {}'.format(lineno))
            if "dblinks" in reaction:
                for i in reaction["dblinks"]:
                    if i[0:2] == "RN":
                        listall = re.split(" ", i[1:])
                        for r in listall:
                            if r[0] == "R":
                                rxn_dict[r] += rxn_mapping[reactions]
    return(rxn_dict)


def gen_transporters(path, gene_asso, tp_sub_dict, tp_fam_dict, chebi_dict,
                     compartment_in, compartment_out):
    with open(os.path.join(path, "transporter_log.tsv"), 'w') \
            as log:
        log.write("compounds for which transporters were predicted, but"
                  " are not found in the model:\n")
        eq = defaultdict(lambda: [])
        for id, genes in gene_asso.items():
            # if there is just one substrate and it is in the porin family
            # assume passive diffusion/porin behavior
            if len(tp_sub_dict[id]) == 1:  # id[0] == "1" and
                for cpd in tp_sub_dict[id]:
                    if cpd in chebi_dict:
                        eq[id].append(("{}_diff".format(chebi_dict[cpd][0]),
                                       Reaction(Direction.Both,
                                       {Compound(chebi_dict[cpd][0]).
                                        in_compartment(compartment_in): -1,
                                        Compound(chebi_dict[cpd][0]).
                                        in_compartment(compartment_out): 1})))
                    else:
                        log.write("{}\t{}\n".format("Compound not in "
                                  "model: ", cpd, id))
            # Otherwise, create a template reaction entry with no
            # formulation for user curation
            else:
                eq[id].append(("TP_{}".format(id),
                              Reaction(Direction.Both,
                                       {Compound("").in_compartment(
                                        compartment_in): -1,
                                        Compound("").in_compartment(
                                        compartment_out): 1})))

        # Sets up the yaml object for the reactions and writes
        # out the parsed reaction information to a reactions.yaml
        # file
        yaml.add_representer(OrderedDict, dict_representer)
        yaml.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
                             dict_constructor)
        yaml_args = {'default_flow_style': False,
                     'encoding': 'utf-8',
                     'allow_unicode': True}

        # construct a new reactions.yaml file based on this.
        reaction_entry_list = []
        mark = FileMark(None, 0, 0)
        rhea = None
        with open(os.path.join(path, "transporters.yaml"),
                  "w") as f:
            for tcdb in eq:
                for rxn in eq[tcdb]:
                    fam_id = ".".join(tcdb.split(".")[0:3])
                    reaction = {'transport_name': rxn[0], 'entry': [rxn[0]],
                                'name': [tp_fam_dict[fam_id]],
                                'equation': [str(rxn[1])],
                                'enzyme': [tcdb],
                                'tcdb_family': tp_fam_dict[fam_id],
                                'substrates': tp_sub_dict[tcdb],
                                'genes': "({})".format(" or ".join(
                                    gene_asso[tcdb]))}
                    entry = ReactionEntry(reaction, rhea, filemark=mark)
                    reaction_entry_list.append(entry)

            yaml.dump(list(model_reactions(reaction_entry_list)), f,
                      **yaml_args)


class main_transporterCommand(Command):
    """Predicts the function of transporter reactions.

    Output based on algorithm outlined in PMID: 18586723"""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument('--annotations', metavar='path',
                            help='Path to the annotation file from eggnog')
        parser.add_argument('--db_substrates', type=str,
                            help='OPTIONAL. Path to a custom substrate file.')
        parser.add_argument('--db_families', type=str,
                            help='OPTIONAL. Path to a custom family file.')
        parser.add_argument('--model', metavar='path',
                            help='Path to the model directory')
        parser.add_argument('--compartment_in', default='c',
                            help='abbreviation for the internal compartment')
        parser.add_argument('--compartment_out', default='e',
                            help='abbreviation for the external compartment')
        parser.add_argument('--col', default=None, help='If providing your own'
                            ' annotation table, specify column for the R, EC, '
                            'or KO number. The default is to parse eggnog'
                            ' output.', type=int)
        super(main_transporterCommand, cls).init_parser(parser)

    def run(self):
        if sys.version_info.minor < 6:
            raise VersionError("Biopython only compatible with python > 3.5.")
        '''Entry point for the transporter assignment'''

        # Check the validity of the input values
        if not self._args.annotations:
            raise InputError('Please specify a path to the annotations')
        if self._args.db_substrates:
            logger.info("Using custom transporter families")
            substrate = self._args.db_substrates
        else:
            logger.info("Using the default transporter families from TCDB. "
                        "Downloaded from: ")
            logger.info("http://www.tcdb.org/cgi-bin/projectv/public/"
                        "getSubstrates.py")
            substrate = resource_filename('psamm',
                                          'external-data/tcdb_substrates.tsv')
        if self._args.db_families:
            logger.info("Using custom transporter families")
            family = self._args.db_families
        else:
            logger.info("Using the default transporter families from TCDB. "
                        "Downloaded from: ")
            logger.info("http://www.tcdb.org/cgi-bin/projectv/public/"
                        "families.py")
            family = resource_filename('psamm',
                                       'external-data/tcdb_families.tsv')
        if not self._args.model:
            raise InputError('Please specify a directory to an existing model')

        # read in the model and build a dictionary of Chebi IDs
        mr = ModelReader.reader_from_path(self._args.model)
        nm = mr.create_model()

        chebi_dict = defaultdict(lambda: [])
        for cpd in nm.compounds:
            if 'ChEBI' in cpd.__dict__['_properties']:
                chebi = re.split(' ', cpd.__dict__['_properties']['ChEBI'])
                for i in chebi:
                    chebi_dict["CHEBI:{}".format(i)].append(cpd.id)

        # Read in the reaction substrates
        tp_sub_dict = defaultdict(lambda: [])
        with open(substrate, 'r', encoding="utf8") as infile:
            for line in infile:
                line = line.rstrip()
                listall = re.split("\t", line)
                substrates = re.split(r"\|", listall[1])
                sub_out = []
                for i in substrates:
                    sub_out.append(re.split(";", i)[0])
                tp_sub_dict[listall[0]] = sub_out

        # read in the reaction families
        tp_fam_dict = defaultdict(lambda: '')
        with open(family, 'r', encoding="utf8") as infile:
            for line in infile:
                line = line.rstrip()
                listall = re.split("\t", line)
                tp_fam_dict[listall[0]] = listall[1]

        # build a gene association dictionary
        # Launch the parse_orthology function to contruct a gene
        # association file
        gene_asso = parse_orthology(self._args.annotations, "tcdb",
                                    self._args.col)

        # Build the transporter databse based on what is in
        # the annotations
        if "model.yaml" in self._args.model:
            path = self._args.model.rstrip("model.yaml")
        else:
            path = self._args.model
        gen_transporters(path, gene_asso, tp_sub_dict, tp_fam_dict,
                         chebi_dict, self._args.compartment_in,
                         self._args.compartment_out)


class main_databaseCommand(Command):
    """Generate a database of compounds and reactions"""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument('--annotations', metavar='path',
                            help='Path to the annotation file from Eggnog')
        parser.add_argument('--type', type=str,
                            help='Define whether to build the model on'
                            ' reaction ID, KO, or EC.\noptions are: '
                            '[R, KO, EC]')
        parser.add_argument('--out', metavar='out',
                            help='Path to the output location for the model '
                            'directory. This will be created for you.')
        parser.add_argument('--col', default=None, help='If providing your '
                            'own annotation table, specify column for the R, '
                            'EC, or KO number. The default is to parse eggnog'
                            ' output.', type=int)
        parser.add_argument('--verbose', help='Report progress  verbose mode',
                            action='store_true')
        parser.add_argument('--default_compartment',
                            help='Define abbreviation for the default '
                            'compartment',
                            default='c')
        parser.add_argument('--rhea', help='Resolve protonation states '
                            'using major microspecies at pH 7.3 using '
                            'Rhea-ChEBI mappings', action='store_true')
        parser.add_argument('--force', help='force rewrite of existing '
                            'directories', action='store_true')
        super(main_databaseCommand, cls).init_parser(parser)

    def run(self):
        """Entry point for the database generation script"""
        # check if required packages are installed
        if 'Bio.KEGG.REST' not in sys.modules:
            quit('No biopython package found. '
                 'Please run <pip install biopython>''')
        if "libchebipy._chebi_entity" not in sys.modules:
            quit('The Chebi API is not installed. '
                 'Please run <pip install libchebipy>')
        if sys.version_info.minor < 6:
            raise VersionError("Biopython only compatible with python > 3.5.")

        # Check the validity of the input values
        if not self._args.annotations:
            raise InputError('Please specify a path to the eggnog annotations')
        if not self._args.type:
            raise InputError('Please specify one of R, KO, or EC as the type')
        if not self._args.out:
            raise InputError('Please specify an output directory')
        if os.path.exists(self._args.out) and not self._args.force:
            raise InputError('The output directory already exists! Exiting...')
            exit()
        elif os.path.exists(self._args.out) is False:
            os.mkdir(self._args.out)

        # Check the format of the eggnog annotation file.
        # Add some code here to check the input file.

        # Launch the parse_orthology function to contruct a
        # gene association file
        ortho_dict = parse_orthology(self._args.annotations, self._args.type,
                                     self._args.col)

        # convert EC to reactions
        if self._args.type == "EC":
            ortho_dict = parse_rxns_from_EC(ortho_dict, self._args.out,
                                            self._args.verbose)
        elif self._args.type == "KO":
            ortho_dict = parse_rxns_from_KO(ortho_dict, self._args.out,
                                            self._args.verbose)

        # Create the model using the kegg api
        create_model_api(self._args.out, ortho_dict, self._args.verbose,
                         self._args.rhea, self._args.default_compartment)


def main(command_class=None, args=None):
    """Run the command line interface with the given :class:`Command`.

    If no command class is specified the user will be able to select a specific
    command through the first command line argument. If the ``args`` are
    provided, these should be a list of strings that will be used instead of
    ``sys.argv[1:]``. This is mostly useful for testing.
    """

    # Set up logging for the command line interface
    if 'PSAMM_DEBUG' in os.environ:
        level = getattr(logging, os.environ['PSAMM_DEBUG'].upper(), None)
        if level is not None:
            logging.basicConfig(level=level)
    else:
        logging.basicConfig(
            level=logging.INFO, format='%(levelname)s: %(message)s')

    parser = argparse.ArgumentParser(description='Options to generate '
                                     'a metabolic model')

    if command_class is not None:
        # Command explicitly given, only allow that command
        command_class.init_parser(parser)
        parser.set_defaults(command=command_class)
    else:
        # Discover all available options
        commands = {}
        for entry in pkg_resources.iter_entry_points(
                'psamm.generate_model'):
            canonical = entry.name.lower()
            if canonical not in commands:
                command_class = entry.load()
                commands[canonical] = command_class
            else:
                logger.warning('Option {} was found more than once!'.format(
                    canonical.name))

        # Create parsers for subcommands
        subparsers = parser.add_subparsers(title='Commands',
                                           metavar='command')
        for name, command_class in sorted(iteritems(commands)):
            title, _, _ = command_class.__doc__.partition('\n\n')
            subparser = subparsers.add_parser(
                name, help=title.rstrip('.'),
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description=_trim(command_class.__doc__))
            subparser.set_defaults(command=command_class)
            command_class.init_parser(subparser)

    parsed_args = parser.parse_args(args)
    command = parsed_args.command(parsed_args)
    try:
        command.run()
    except CommandError as e:
        parser.error(text_type(e))


class ReactionEntry(object):
    """Representation of entry in KEGG reaction file"""

    def __init__(self, values, rhea, filemark=None):
        self.values = dict(values)
        if 'entry' not in values:
            raise ParseError2('Missing reaction identifier')
        if 'transport_name' in values:
            self._id = values['transport_name']
        else:
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
    def tcdb_family(self):
        if 'tcdb_family' not in self.values:
            return None
        return self.values['tcdb_family']

    @property
    def substrates(self):
        if 'substrates' not in self.values:
            return None
        return self.values['substrates']

    @property
    def genes(self):
        if 'genes' not in self.values:
            return None
        return self.values['genes']

    @property
    def orthology(self):
        if 'orthology' in self.values:
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


class CompoundEntry(object):
    """Representation of entry in KEGG compound file"""

    def __init__(self, values, rhea, filemark=None):
        self.values = dict(values)
        if 'entry' not in values:
            raise ParseError2('Missing compound identifier')
        self._id, _ = values['entry'][0].split(None, 1)
        self._filemark = filemark
        self._charge = None
        self._chebi = None
        self._chebi_all = None
        self.rhea = rhea
        self.initialize_charge()

    def initialize_charge(self):
        """
        Sets the _charge, _chebi, and _chebi_all attributes
        'rhea_db' is initialized as a global in generate_model_api
        if --rhea is supplied this funcion looks for rhea_db in the
        global namespace decide if rhea is used
        --- Logic for selecting the best chebi ID ---
        if not using rhea:
            use the first chebi ID given by KEGG
        elif using rhea:
            if all KEGG-chebi IDs map to same ID in rhea:
                use the single ID
            elif KEGG-chebi IDs map to different IDs in rhea:
                use the first chebi ID given by KEGG
            elif the KEGG-chebi IDs don't have mappings in rhea:
                use the first chebi ID given by KEGG
        """
        if self.rhea is not None:
            use_rhea = True
            rhea_db = self.rhea
        else:
            use_rhea = False
        for DB, ID in self.dblinks:
            if DB == "ChEBI":
                id_list = ID.split(" ")
                if use_rhea:
                    rhea_id_list = rhea_db.select_chebi_id(id_list)
                    if len(rhea_id_list) == 0:  # no chebis map to rhea
                        self._chebi = id_list[0]
                        self._chebi_all = id_list
                    elif len(set(rhea_id_list)) == 1:  # chebi map to same rhea
                        self._chebi = rhea_id_list[0]
                        self._chebi_all = set(id_list + [rhea_id_list[0]])
                    else:  # chebis map to different rheas
                        self._chebi = id_list[0]
                        self._chebi_all = set(id_list + rhea_id_list)
                else:  # --rhea not given
                    self._chebi = id_list[0]
                    self._chebi_all = list(id_list)

        # libchebipy update charge and formula
        if self._chebi is not None:
            this_chebi_entity = ChebiEntity(self._chebi)
            try:
                try:
                    # libchebipy sometimes fails with an index error
                    # on the first time running. We have not been able
                    # to fix the source of this error, but catching the
                    # index error and repeating the operation appears to
                    # fix this
                    self._charge = int(this_chebi_entity.get_charge())
                    self.values['formula'] = [this_chebi_entity.get_formula()]
                except IndexError:
                    self._charge = int(this_chebi_entity.get_charge())
                    self.values['formula'] = [this_chebi_entity.get_formula()]
            except ValueError:  # chebi entry has no charge; leave as None
                pass

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
    def dblinks(self):
        if 'dblinks' in self.values:
            for line in self.values['dblinks']:
                database, entry = line.split(':', 1)
                yield database.strip(), entry.strip()
        else:
            return None

    @property
    def charge(self):
        return self._charge

    @property
    def chebi(self):
        return self._chebi

    @property
    def chebi_all(self):
        if self._chebi_all is not None:
            return ', '.join(self._chebi_all)
        else:
            return None

    @property
    def comment(self):
        if 'comment' not in self.values:
            return None
        try:
            return '\n'.join(self.values['comment'])
        except TypeError:
            return self.values['comment']

    def __getitem__(self, name):
        if name not in self.values:
            raise AttributeError('Attribute does not exist: {}'.format(name))
        return self._values[name]

    @property
    def filemark(self):
        return self._filemark

    def __repr__(self):
        return '<CompoundEntry "{}">'.format(self.id)


class RheaDb(object):
    """Allows storing and searching Rhea db"""

    def __init__(self, filepath):
        self._values = self._parse_db_from_tsv(filepath)

    def _parse_db_from_tsv(self, filepath):
        db = {}
        with open(filepath, 'r') as f:
            for line in f:
                split = line.split('\t')
                db[split[0]] = split[1]
        return db

    def select_chebi_id(self, id_list):
        rhea_id_list = []
        for x in id_list:
            try:
                rhea_id_list.append(self._values[x])
            except KeyError:
                pass
        return rhea_id_list
