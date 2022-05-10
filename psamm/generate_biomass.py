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
# Copyright 2021-2022  Jason Vailionis <jason_vailionis@uri.edu>

from psamm.datasource.reaction import Reaction, Compound, Direction
from psamm.datasource import native
from collections import Counter
from six import add_metaclass
from pkg_resources import resource_filename
import logging
import sys
import os
import abc
import yaml
import pandas as pd
import numpy as np
if sys.version_info.minor > 5:
    try:
        from Bio import SeqIO
        from Bio import Seq
    except ImportError:
        quit("No biopython package found. "
             "Please run <pip install biopython>")

logger = logging.getLogger(__name__)


@add_metaclass(abc.ABCMeta)
class Command(object):
    """Represents a command in the interface, operating on a model.

    The constructor will be given the NativeModel and the command line
    namespace. The subclass must implement :meth:`run` to handle command
    execution. The doc string will be used as documentation for the command
    in the command line interface.

    In addition, :meth:`init_parser` can be implemented as a classmethod which
    will allow the command to initialize an instance of
    :class:`argparse.ArgumentParser` as desired. The resulting argument
    namespace will be passed to the constructor.
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
            logger.debug("Command failure caused by exception!", exc_info=exc)
        sys.exit(1)


class InputError(Exception):
    """Exception used to signal a general input error."""


class CommandError(Exception):
    """Error from running a command.

    This should be raised from a ``Command.run()`` if any arguments are
    misspecified. When the command is run and the ``CommandError`` is raised,
    the caller will exit with an error code and print appropriate usage
    information.
    """


# Returns df with additional columns for stoichiometry
def calc_stoichiometry_df(df, counts):
    df["counts"] = pd.Series(counts)
    df["composition"] = df.counts / df.counts.sum()
    df["mass"] = df.mol_weight*df.composition
    df["stoichiometry"] = df.composition / df.mass.sum()
    df = df.fillna(0)
    return df


def load_compound_data():
    df = pd.read_csv(resource_filename("psamm",
                     "external-data/biomass_compound_descriptions.tsv"),
                     sep="\t", na_values=["None"])
    df.set_index(df.id.values, inplace=True)
    return df


def make_template(compounds_df):
    template = compounds_df[["id", "name"]]
    template["custom_id"] = np.nan
    return template


def apply_config(df, config_path):
    config = pd.read_csv(config_path)
    save_index = df.index
    df = pd.merge(df, config, how="left", on=["id", "name"])
    df.set_index(save_index, inplace=True)
    df.custom_id = df.custom_id.fillna(df.id)
    df.id = df.custom_id
    return df

# return model_dict


def load_model(model_path):
    with open(model_path, "r") as f:
        model_dict = yaml.safe_load(f)
    return model_dict


def check_missing_cpds(model_path, df, config, yaml_args):
    # Check if needed compounds are in the model, add them if not
    model_dir = os.path.dirname(model_path)
    model_reader = native.ModelReader.reader_from_path(model_path)
    compounds = [r.id for r in model_reader.parse_compounds()]
    if len(compounds) > 0:
        missing_cpds = df.query("id not in @compounds")
        droplist = ["type", "code"]
        if config:
            droplist.append("custom_id")
        missing_cpds = missing_cpds.drop(columns=droplist)
        if len(missing_cpds) > 0:
            logger.warning("\033[1mThe following compounds are required for "
                           "biomass equations but are missing from the model: "
                           "\033[0m{}\n\033[1mPredefined versions of these "
                           "compounds will be generated in {}. To use custom "
                           "biomass compound IDs, use --config\033[0m"
                           "".format(list(missing_cpds.id),
                                     os.path.join(model_dir,
                                                  "biomass_compounds.yaml")))
            missing_cpd_entries = [v.dropna().to_dict()
                                   for k, v in missing_cpds.iterrows()]
            with open(os.path.join(model_dir,
                                   "biomass_compounds.yaml"), "a+") as f:
                yaml.dump(missing_cpd_entries, f, **yaml_args)


def fix_cell_compartment(model_dict):
    try:
        cell_compartment = model_dict["default_compartment"]
    except KeyError:
        cell_compartment = "c"
        model_dict["default_compartment"] = "c"
    if cell_compartment is None:
        cell_compartment = "c"
        model_dict["default_compartment"] = "c"
    return cell_compartment


def update_model_yaml(model_path, model_dict, cell_compartment,
                      biomass_name, yaml_args):
    exists = False  # Checking if biomass_reactions.yaml already included
    for file in model_dict["reactions"]:
        if "biomass_reactions.yaml" in file["include"]:
            exists = True
    if not exists:
        model_dict["reactions"].append({"include":
                                        "./biomass_reactions.yaml"})

    exists = False  # Checking if biomass_compounds.yaml already included
    for file in model_dict["compounds"]:
        if "biomass_compounds.yaml" in file["include"]:
            exists = True
    if not exists:
        model_dict["compounds"].append({"include":
                                        "./biomass_compounds.yaml"})

    model_dict["biomass"] = biomass_name

    os.rename(model_path, model_path + ".tmp")
    with open(model_path, "w") as f:
        yaml.dump(model_dict, f, **yaml_args)
    os.remove(model_path + ".tmp")


def count_DNA(genome, df, cell_compartment, decimals=6):
    DNA = df[df.type == "DNA"].copy().set_index("code")
    DNA_counts = Counter()
    for seq in genome:
        DNA_counts += Counter(genome[seq].seq)

    calc_stoichiometry_df(DNA, DNA_counts)
    total = DNA.stoichiometry.sum()
    ids_fwd = list(DNA.id) + list(df.loc[["C00002", "C00001"]].id.values)
    stoich_fwd = list(DNA.stoichiometry) + [2*total, 2*total]
    compound_fwd = {Compound(id).in_compartment(cell_compartment):
                    round(-1*n, decimals) for id, n in zip(ids_fwd,
                                                           stoich_fwd)}
    ids_rev = list(df.loc[["dna", "C00008", "C00009", "C00013"]].id.values)
    stoich_rev = [1, 2*total, 2*total, total]
    compound_rev = {Compound(id).in_compartment(cell_compartment):
                    round(n, decimals) for id, n in zip(ids_rev, stoich_rev)}
    dna_rxn = Reaction(Direction.Forward, {**compound_fwd, **compound_rev})
    dna_rxn_entry = {"id": "dna_met",
                     "name": "DNA",
                     "equation": str(dna_rxn),
                     "pathways": ["Biomass"]}
    return dna_rxn_entry


def count_Prot(proteome, df, cell_compartment, decimals=6):
    Prot = df[df.type == "Prot"].copy().set_index("code")
    Prot_counts = Counter()
    for seq in proteome:
        Prot_counts += Counter(proteome[seq].seq)

    # Protein Reaction formation
    calc_stoichiometry_df(Prot, Prot_counts)
    total = Prot.stoichiometry.sum()
    ids_fwd = list(df[df.type == "aatrna"].id) + \
        list(df.loc[["C00044", "C00001"]].id.values)
    stoich_fwd = list(Prot.stoichiometry) + [2*total, 2*total]
    compound_fwd = {Compound(id).in_compartment(cell_compartment):
                    round(-1*n, decimals) for id, n in zip(ids_fwd,
                                                           stoich_fwd)}
    ids_rev = list(df[df.type == "trna"].id) +\
        list(df.loc[["protein", "C00035", "C00009"]].id.values)
    stoich_rev = list(Prot.stoichiometry) + [1, 2*total, 2*total]
    compound_rev = {Compound(id).in_compartment(cell_compartment):
                    round(n, decimals) for id, n in zip(ids_rev, stoich_rev)}
    prot_rxn = Reaction(Direction.Forward, {**compound_fwd, **compound_rev})
    prot_rxn_entry = {"id": "pro_met",
                      "name": "Protein",
                      "equation": str(prot_rxn),
                      "pathways": ["Biomass"]}
    return prot_rxn_entry


def count_RNA(genome, gff, df, cell_compartment, decimals=6):
    RNA = df[df.type == "RNA"].copy().set_index("code")
    RNA_counts = Counter()
    # Parsing through the gff and counting bases in coding regions to get
    # the total RNA base counts
    with open(gff, "r") as gff_file:
        for line in gff_file:
            if not line.startswith("#"):
                L = line.strip().split("\t")
                if L[2] == "CDS":
                    try:
                        seq_entry = genome[L[0]]
                    except KeyError:
                        logger.warning("The annotation {s} does not match "
                                       "any sequence names in genome file, "
                                       "ignoring...".format(s=L[0]))
                    # FIX STRANDEDNESS: if - use reverse complement
                    if L[6] == "-":
                        RNA_counts += Counter(Seq.complement(
                            seq_entry.seq[int(L[3]):int(L[4])]))
                    else:
                        RNA_counts += Counter(
                            seq_entry.seq[int(L[3]):int(L[4])])
    try:
        RNA_counts["U"] = RNA_counts.pop("T")  # renames T to U for RNA
    except KeyError:
        pass  # If there are no T's, pass

    # RNA Reaction formation
    calc_stoichiometry_df(RNA, RNA_counts)
    total = RNA.stoichiometry.sum()
    ids_fwd = list(RNA.id) + list(df.loc[["C00001"]].id.values)
    stoich_fwd = [RNA.loc["A"].stoichiometry + 2*total] +\
        list(RNA.stoichiometry)[1:] + [2*total]
    compound_fwd = {Compound(id).in_compartment(cell_compartment):
                    round(-1*n, decimals) for id, n in zip(ids_fwd,
                                                           stoich_fwd)}
    ids_rev = list(df.loc[["rna", "C00008", "C00009", "C00013"]].id.values)
    stoich_rev = [1, 2*total, 2*total, total]
    compound_rev = {Compound(id).in_compartment(cell_compartment):
                    round(n, decimals) for id, n in zip(ids_rev, stoich_rev)}
    rna_rxn = Reaction(Direction.Forward, {**compound_fwd, **compound_rev})
    rna_rxn_entry = {"id": "rna_met",
                     "name": "RNA",
                     "equation": str(rna_rxn),
                     "pathways": ["Biomass"]}
    return rna_rxn_entry


def return_biomass_rxn(df, biomass_name, cell_compartment):
    # Biomass reaction formation
    compound_fwd = {Compound(id).in_compartment(cell_compartment): -1
                    for id in list(df.loc[["dna",
                                           "rna", "protein"]].id.values)}
    compound_rev = {
        Compound(df.loc["biomass"].id).in_compartment(cell_compartment): 1}
    bio_rxn = Reaction(Direction.Forward, {**compound_fwd, **compound_rev})

    bio_rxn_entry = {"id": biomass_name,
                     "name": "Biomass",
                     "equation": str(bio_rxn),
                     "pathways": ["Biomass"]}
    return bio_rxn_entry


def return_bio_sink_rxn(df, cell_compartment):
    # Biomass sink reaction
    bio_sink = Reaction(Direction.Forward, {
        Compound(df.loc["biomass"].id).in_compartment(cell_compartment): -1})

    bio_sink_entry = {"id": "sink_biomass",
                      "name": "Biomass accumulation",
                      "equation": str(bio_sink),
                      "pathways": ["Biomass"]}
    return bio_sink_entry


def return_trna_rxns(df, cell_compartment):
    # tRNA charging reactions
    tRNA_rxns = []
    df_reac = pd.read_csv(
        resource_filename("psamm",
                          "external-data/biomass_reaction_descriptions.tsv"),
        sep="\t")
    for index, row in df_reac.iterrows():
        fwd, rev = row.equation.split(" <=> ")
        ids_fwd = fwd.split(" + ")
        ids_rev = rev.split(" + ")
        compound_fwd = {Compound(id).in_compartment(cell_compartment): -1
                        for id in list(df.loc[ids_fwd].id)}
        compound_rev = {Compound(id).in_compartment(cell_compartment): 1
                        for id in list(df.loc[ids_rev].id)}
        tRNA_rxn = Reaction(Direction.Both, {**compound_fwd, **compound_rev})
        tRNA_rxn_entry = {"id": row.id,
                          "name": row["name"],
                          "equation": str(tRNA_rxn),
                          "enzyme": row["enzymes"],
                          "pathways": ["Biomass"]}
        tRNA_rxns.append(tRNA_rxn_entry)
    return tRNA_rxns


class main(Command):
    """Generate a database of compounds and reactions"""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument("--genome", metavar="path",
                            help="[REQUIRED] Path to the genome "
                            "in fasta format")
        parser.add_argument("--proteome", metavar="path",
                            help="[REQUIRED] Path to the proteome "
                            "in fasta format")
        parser.add_argument("--gff", metavar="path",
                            help="[REQUIRED] Specify path to gff containing "
                            "transcript annotations. Only annotations "
                            "specified as \"CDS\" in the third column will be "
                            "used. Annotations must correspond to sequences in"
                            " the --genome file.")
        parser.add_argument("--model", metavar="path",
                            help="[REQUIRED] Path to the model file")
        parser.add_argument("--biomass", metavar="name", default="biomass",
                            help="Name to use for the biomass reaction, "
                            "default 'biomass'")
        parser.add_argument("--config", metavar="path",
                            help="If you are using custom compound names for "
                            "the biomass-related compounds (e.g. amino acids, "
                            "DNA bases), specify them via a 3-column config "
                            "file. Use --generate-config to create a template "
                            "file.")
        parser.add_argument("--generate-config", action="store_true",
                            help="Makes a template config file for --config. "
                            "Compound IDs added in the third column will "
                            "replace the default kegg IDs for the compound.")
        super(main, cls).init_parser(parser)

    def run(self):
        """Entry point for the biomass reaction generation script"""

        if not self._args.generate_config:
            if not self._args.genome:
                raise InputError("Please specify a path to the genome "
                                 "in fasta format with --genome")
            if not self._args.proteome:
                raise InputError("Specify path to proteome in fasta format "
                                 "with --proteome")
            if not self._args.gff:
                raise InputError("Specify path to gff containing transcript "
                                 "annotations with --gff.")
            if not self._args.model:
                raise InputError("Please specify the path to an existing "
                                 "model with --model")
            if os.path.isdir(self._args.model):
                if os.path.exists(os.path.join(self._args.model,
                                               "model.yaml")):
                    model_path = os.path.join(self._args.model, "model.yaml")
                else:
                    raise InputError("No model.yaml file found!")
            else:
                model_path = self._args.model

        # SETTING UP DATA
        pd.options.mode.chained_assignment = None
        df = load_compound_data()
        if self._args.generate_config:
            template = make_template(df)
            template.to_csv(sys.stdout, index=False)
            quit()

        faa = self._args.proteome
        fna = self._args.genome
        gff = self._args.gff
        model_dir = os.path.dirname(model_path)

        yaml_args = {"default_flow_style": False,
                     "sort_keys": False,
                     "encoding": "utf-8",
                     "allow_unicode": True,
                     "width": float("inf")}

        # Handling the config file
        if self._args.config:
            df = apply_config(df, self._args.config)

        # Adds missing compounds, updates model.yaml
        # Makes biomass compounds and reactions yamls
        logger.info("Checking for existing compounds in the model")
        check_missing_cpds(model_path, df, self._args.config, yaml_args)
        model_dict = load_model(model_path)
        cell_compartment = fix_cell_compartment(model_dict)

        logger.info("Updating {}".format(os.path.abspath(model_path)))
        update_model_yaml(model_path, model_dict, cell_compartment,
                          self._args.biomass, yaml_args)

        # Counting bases/amino acids and making reaction entries
        genome = {seq.name: seq.upper() for seq in SeqIO.parse(fna, "fasta")}
        proteome = {seq.name: seq.upper() for seq in SeqIO.parse(faa, "fasta")}

        logger.info("Calculating stoichiometry for DNA, RNA, and amino acids")
        dna_rxn_entry = count_DNA(genome, df, cell_compartment)
        prot_rxn_entry = count_Prot(proteome, df, cell_compartment)
        rna_rxn_entry = count_RNA(genome, gff, df, cell_compartment)

        # Making the reaction entries for other biomass rxns
        bio_rxn_entry = return_biomass_rxn(
            df, self._args.biomass, cell_compartment)
        bio_sink_entry = return_bio_sink_rxn(df, cell_compartment)
        tRNA_rxns = return_trna_rxns(df, cell_compartment)

        # GENERATING biomass_reactions.yaml
        logger.info("Generating new biomass reactions in "
                    "{}/biomass_reactions.yaml".format(
                        os.path.abspath(model_dir)))

        with open(os.path.join(model_dir, "biomass_reactions.yaml"), "w") as f:
            yaml.dump([dna_rxn_entry, rna_rxn_entry,
                       prot_rxn_entry, bio_rxn_entry,
                       bio_sink_entry] + tRNA_rxns,
                      f, **yaml_args)
