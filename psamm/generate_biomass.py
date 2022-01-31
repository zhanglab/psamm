from psamm.datasource.reaction import Reaction, Compound, Direction
from psamm.datasource import native
from collections import Counter
from collections import OrderedDict
from six import add_metaclass, iteritems, itervalues, text_type, PY3
from pkg_resources import resource_filename
import logging
import sys
import os
import abc
import yaml
import pandas as pd
import numpy as np
try:
    from Bio import SeqIO
    from Bio import Seq
except ModuleNotFoundError:
    quit("No biopython package found. Please run <pip install biopython>")

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


class main(Command):
    """Generate a database of compounds and reactions"""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument("--genome", metavar="path",
            help="[REQUIRED] Path to the genome in fasta format")
        parser.add_argument("--proteome", metavar="path",
            help="[REQUIRED] Path to the proteome in fasta format")
        parser.add_argument("--gff", metavar = "path",
            help="[REQUIRED] Specify path to gff containing transcript \
                  annotations. Only annotations specified as \"CDS\" in the \
                  third column will be used. Annotations must correspond to \
                  sequences in the --genome file.")
        parser.add_argument("--model", metavar="path",
            help="[REQUIRED] Path to the model file")
        parser.add_argument("--biomass", metavar="name", default = "biomass",
            help="Name to use for the biomass reaction, default 'biomass'")
        parser.add_argument("--config", metavar="path",
            help="If you are using custom compound names for the biomass-"
                 "related compounds (e.g. amino acids, DNA bases), specify "
                 "them via a 3-column config file. Use --generate-config "
                 "to create a template file.")
        parser.add_argument("--generate-config", action = "store_true",
            help="Makes a template config file for --config. Compound IDs "
                 "added in the third column will replace the default kegg IDs "
                 "for the compound")
        parser.add_argument("--non-bacteria", action = "store_true",
            help="Exclude the bacteria-specific formylated methionine reaction")
        ##------------------------ To do---------------------------------------
        ## need to check first if reactions are already in model?
        ## why are the formulas of the trnas 'R' in Cbes; how should I save them?
        ## Make the biomass equations work again
        ## Might need to catch spaces in the config compound names, havent tested
        ## What are R04212 and R03905?
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
                    "annotations with --gff. Only annotations specified as \"CDS\" "
                    "in the third column will be used. Annotations must correspond "
                    "to sequences in the --genome file.")
            if not self._args.model:
                raise InputError("Please specify the path to an existing "
                                 "model with --model")

        ### SETTING UP DATA AND FUNCTIONS ###
        pd.options.mode.chained_assignment = None
        # Returns df with additional columns for stoichiometry
        def calc_stoichiometry_df(df, counts):
            df["counts"] = pd.Series(counts)
            df["composition"] = df.counts / df.counts.sum()
            df["mass"] = df.mol_weight*df.composition
            df["stoichiometry"] = df.composition / df.mass.sum()
            df = df.fillna(0)
            return df

        df = pd.read_csv(resource_filename("psamm",
                            "external-data/biomass_compound_descriptions.tsv"),
                        sep = "\t", na_values = ["None"])
        df.set_index(df.id.values, inplace = True)

        if self._args.generate_config:
            template = df[["id", "name"]]
            template["custom_id"] = np.nan
            template.to_csv(sys.stdout, index = False)
            quit()

        yaml_args = {"default_flow_style": False,
                     "sort_keys": False,
                     "encoding": "utf-8",
                     "allow_unicode": True,
                     "width": float("inf")}

        faa = self._args.proteome
        fna = self._args.genome
        gff = self._args.gff
        model_path = self._args.model
        model_dir = os.path.dirname(model_path)

        ### Handling the config file ###
        if self._args.config:
            config = pd.read_csv(self._args.config)
            save_index = df.index
            df = pd.merge(df, config, how = "left", on = ["id", "name"])
            df.set_index(save_index, inplace = True)
            df.custom_id = df.custom_id.fillna(df.id)
            df.id = df.custom_id

        # Slice the dataframe after config so that new ids are maintained
        Prot = df[df.type == "Prot"].copy().set_index("code")
        DNA = df[df.type == "DNA"].copy().set_index("code")
        RNA = df[df.type == "RNA"].copy().set_index("code")
        #Other = df[df.type == "Other"]

        ### Check if needed compounds are in the model, add them if not
        model_reader = native.ModelReader.reader_from_path(model_path)
        compounds = [r.id for r in model_reader.parse_compounds()]
        missing_cpds = df.query("id not in @compounds")
        droplist = ["type"]
        if self._args.config:
            droplist.append("custom_id")
        missing_cpds = missing_cpds.drop(columns = droplist)
        if len(missing_cpds) > 0:
            logger.warning("\033[1mThe following compounds are required for "
                "biomass equations but are missing from the model: "
                "\033[0m{}\n\033[1mEither create entries for these compounds "
                "or specify corresponding names to existing "
                "compounds using --config. \nPredefined versions of these "
                "compounds will be generated in "
                "{}\033[0m".format(list(missing_cpds.id),os.path.join(model_dir,
                                                    "biomass_compounds.yaml")))
            missing_cpd_entries = [v.dropna().to_dict()
                                    for k,v in missing_cpds.iterrows()]
            with open(os.path.join(model_dir, "biomass_compounds.yaml"), "a+") as f:
                yaml.dump(missing_cpd_entries, f, **yaml_args)

        ### Updating model.yaml ###
        logger.info("Updating model in {}".format(os.path.abspath(model_path)))

        with open(model_path, "r") as f:
            model_dict = yaml.safe_load(f)

        exists = False # Checking if biomass_reactions.yaml already included
        for file in model_dict["reactions"]:
            if "biomass_reactions.yaml" in file["include"]:
                exists = True
        if not exists:
            model_dict["reactions"].append({"include":
                                            "./biomass_reactions.yaml"})

        exists = False # Checking if biomass_compounds.yaml already included
        for file in model_dict["compounds"]:
            if "biomass_compounds.yaml" in file["include"]:
                exists = True
        if not exists:
            model_dict["compounds"].append({"include":
                                            "./biomass_compounds.yaml"})

        model_dict["biomass"] = self._args.biomass
        cell_compartment = model_dict["default_compartment"]
        if cell_compartment is None:
            cell_compartment = "c"
            model_dict["default_compartment"] = "c"

        os.rename(model_path, model_path + ".tmp")
        with open(model_path, "w") as f:
            yaml.dump(model_dict, f, **yaml_args)
        os.remove(model_path + ".tmp")

        ### DATA PROCESSING ###
        logger.info("Counting nucleotides and amino acids")
        DNA_counts = Counter()
        RNA_counts = Counter()
        Prot_counts = Counter()

        # genome and proteome are dicts with "seqname": seqIO_entry
        # Counts look like {"A": 200, "C": 380..}
        genome = {seq.name: seq.upper() for seq in SeqIO.parse(fna, "fasta")}
        proteome = {seq.name: seq.upper() for seq in SeqIO.parse(faa, "fasta")}

        for seq in genome:
            DNA_counts += Counter(genome[seq].seq)

        for seq in proteome:
            Prot_counts += Counter(proteome[seq].seq)

        # Parsing through the gff and counting bases in coding regions to get
        # the total RNA base counts
        with open(gff, "r") as gff_file:
            for line in gff_file:
                if not line.startswith("#"):
                    l = line.strip().split("\t")
                    if l[2] == "CDS":
                        try:
                            seq_entry = genome[l[0]]
                        except KeyError:
                            logger.warning("The annotation {s} does not match "
                                "any sequence names in {f}, ignoring...".format(
                                s = l[0], f = fna
                                ))
                        ## FIX STRANDEDNESS: if - use reverse complement
                        if l[6] == "-":
                            RNA_counts += Counter(Seq.complement(
                                            seq_entry.seq[int(l[3]):int(l[4])]))
                        else:
                            RNA_counts += Counter(
                                            seq_entry.seq[int(l[3]):int(l[4])])

        RNA_counts["U"] = RNA_counts.pop("T") #renames T to U for RNA

        ### DNA Reaction formation ###
        logger.info("Calculating stoichiometry for biomass equations")

        decimals = 6

        calc_stoichiometry_df(DNA, DNA_counts)
        total = DNA.stoichiometry.sum()
        ids_fwd = list(DNA.id) + list(df.loc[["C00002", "C00001"]].id.values)
        stoich_fwd = list(DNA.stoichiometry) + [2*total, 2*total]
        compound_fwd = {Compound(id).in_compartment(cell_compartment):
                    round(-1*n, decimals) for id,n in zip(ids_fwd, stoich_fwd)}
        ids_rev = list(df.loc[["dna", "C00008", "C00009", "C00013"]].id.values)
        stoich_rev = [1, 2*total, 2*total, total]
        compound_rev = {Compound(id).in_compartment(cell_compartment):
                    round(n, decimals) for id,n in zip(ids_rev, stoich_rev)}
        dna_rxn = Reaction(Direction.Forward, {**compound_fwd, **compound_rev})
        dna_rxn_entry = {"id": "dna_met",
                         "name": "DNA",
                         "equation": str(dna_rxn),
                         "pathways": ["Biomass"]}

        ### RNA Reaction formation
        calc_stoichiometry_df(RNA, RNA_counts)
        total = RNA.stoichiometry.sum()
        ids_fwd = list(RNA.id) + list(df.loc[["C00001"]].id.values)
        stoich_fwd = [RNA.loc["A"].stoichiometry + 2*total] +\
                    list(RNA.stoichiometry)[1:] + [2*total]
        compound_fwd = {Compound(id).in_compartment(cell_compartment):
                    round(-1*n, decimals) for id,n in zip(ids_fwd, stoich_fwd)}
        ids_rev = list(df.loc[["rna", "C00008", "C00009", "C00013"]].id.values)
        stoich_rev = [1, 2*total, 2*total, total]
        compound_rev = {Compound(id).in_compartment(cell_compartment):
                        round(n, decimals) for id,n in zip(ids_rev, stoich_rev)}
        rna_rxn = Reaction(Direction.Forward, {**compound_fwd, **compound_rev})
        rna_rxn_entry = {"id": "rna_met",
                         "name": "RNA",
                         "equation": str(rna_rxn),
                         "pathways": ["Biomass"]}

        ### Protein Reaction formation
        calc_stoichiometry_df(Prot, Prot_counts)
        total = Prot.stoichiometry.sum()
        ids_fwd = list(df[df.type == "aatrna"].id) + \
                    list(df.loc[["C00044", "C00001"]].id.values)
        stoich_fwd = list(Prot.stoichiometry) + [2*total, 2*total]
        compound_fwd = {Compound(id).in_compartment(cell_compartment):
                    round(-1*n, decimals) for id,n in zip(ids_fwd, stoich_fwd)}
        ids_rev = list(df[df.type == "trna"].id) +\
                    list(df.loc[["protein", "C00035", "C00009"]].id.values)
        stoich_rev = list(Prot.stoichiometry) + [1, 2*total, 2*total]
        compound_rev = {Compound(id).in_compartment(cell_compartment):
                    round(n, decimals) for id,n in zip(ids_rev, stoich_rev)}
        prot_rxn = Reaction(Direction.Forward, {**compound_fwd, **compound_rev})
        prot_rxn_entry = {"id": "pro_met",
                          "name": "Protein",
                          "equation": str(prot_rxn),
                          "pathways": ["Biomass"]}

        ### Biomass reaction formation ###
        compound_fwd = {Compound(id).in_compartment(cell_compartment): -1
                        for id in list(df.loc[["dna", "rna", "protein"]].id.values)}
        compound_rev = {Compound(df.loc["biomass"].id).in_compartment(cell_compartment): 1}
        bio_rxn = Reaction(Direction.Forward, {**compound_fwd, **compound_rev})

        bio_rxn_entry = {"id": self._args.biomass,
                         "name": "Biomass",
                         "equation": str(bio_rxn),
                         "pathways": ["Biomass"]}

        ### Biomass sink reaction ###
        bio_sink = Reaction(Direction.Forward, {
                    Compound(df.loc["biomass"].id).in_compartment(cell_compartment): -1})

        bio_sink_entry = {"id": "sink_biomass",
                         "name": "Biomass accumulation",
                         "equation": str(bio_sink),
                         "pathways": ["Biomass"]}

        #### tRNA charging reactions ###
        tRNA_rxns = []
        df_reac = pd.read_csv(resource_filename("psamm",
                            "external-data/biomass_reaction_descriptions.tsv"),
                            sep = "\t")
        if self._args.non_bacteria:
            df_reac = df_reac[df_reac["id"].str.contains("R03940") == False]
        for index, row in df_reac.iterrows():
            fwd, rev = row.equation.split(" <=> ")
            ids_fwd = fwd.split(" + ")
            ids_rev = rev.split(" + ")
            #list(df.loc[ids_fwd].id)
            compound_fwd = {Compound(id).in_compartment(cell_compartment): -1
                            for id in list(df.loc[ids_fwd].id)}
            compound_rev = {Compound(id).in_compartment(cell_compartment): 1
                            for id in list(df.loc[ids_rev].id)}
            tRNA_rxn = Reaction(Direction.Both, {**compound_fwd,**compound_rev})
            tRNA_rxn_entry = {"id": row.id,
                              "name": row["name"],
                              "equation": str(tRNA_rxn),
                              "enzyme": row["enzymes"],
                              "pathways": ["Biomass"]}
            tRNA_rxns.append(tRNA_rxn_entry)

        ### GENERATING biomass_reactions.yaml ###
        logger.info("Generating new biomass reactions in "
                "{}/biomass_reactions.yaml".format(os.path.abspath(model_dir)))

        with open(os.path.join(model_dir, "biomass_reactions.yaml"), "w") as f:
            yaml.dump([dna_rxn_entry, rna_rxn_entry,
                       prot_rxn_entry, bio_rxn_entry,
                       bio_sink_entry] + tRNA_rxns,
                       f, **yaml_args)
