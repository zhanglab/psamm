from psamm.datasource.reaction import Reaction, Compound, Direction
from psamm.datasource import native
from collections import Counter
from collections import OrderedDict
from six import add_metaclass, iteritems, itervalues, text_type, PY3
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
            help="Path to the genome in fasta format")
        parser.add_argument("--proteome", metavar="path",
            help="Path to the proteome in fasta format")
        parser.add_argument("--gff", metavar = "path",
            help="Specify path to gff containing transcript \
                  annotations. Only annotations specified as \"CDS\" in the \
                  third column will be used. Annotations must correspond to \
                  sequences in the --genome file.")
        parser.add_argument("--model", metavar="path",
            help="Path to the model directory")
        parser.add_argument("--biomass", metavar="name", default = "biomass",
            help="Name of the biomass reaction, default 'biomass'")
        parser.add_argument("--config", metavar="path",
            help="If you are using custom compound names for the biomass-"
                 "related compounds (e.g. amino acids, DNA bases), specify "
                 "them via a 3-column config file. Use --generate-config "
                 "to create a template file.")
        parser.add_argument("--generate-config", action = "store_true",
            help="Makes a template config file for --config. Compound IDs "
                 "added in the third column will replace the default kegg IDs "
                 "for the compound")
        ##------------------------ To do---------------------------------------
        ## Add tRNA charging reactions
        ##  -> use download_kegg_api to get the reaction info or store it?
        ##  -> need to check first if those reactions are already in model
        ##  -> add formylated methionine reaction too (R03940)
        ## Check if required compounds are in the model currently
        ##  -> add compounds to compounds_biomass.yaml? store kegg information or download like reactions?
        ## Add new compounds to new compounds file
        ## Check if dna/rna/protein/biomass are already in compounds.yaml before adding
        ### Combine the dataframes into one with index being the kegg id
        ### and a second index being DNA/RNA/Prot/other
        ### so that the config can be generated without concatening tables
        ### and so that the edited config can be joined in a single merge.
        ### save slices of the dataframe as DNA/RNA/Prot/other so no refactoring
        ### There will be another table for the tRNA charge reactions
        ### which relates the 3 sets of IDs and will allow me to include the
        ### formylated methionine reaction and cpds more easily
        ## Read dataframe from file instead of writing it into the code
        ## So many of these charges disagree with chebi..
        ## e.g. C00131 which is -4 in Cbes model but 0 on chebi
        ## why are the formulas of the trnas 'R' in Cbes; how should I save them?
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
                raise InputError("Specify path to gff containing tggranscript "
                    "annotations with --gff. Only annotations specified as \"CDS\" "
                    "in the third column will be used. Annotations must correspond "
                    "to sequences in the --genome file.")
            if not self._args.model:
                raise InputError("Please specify the path to an existing "
                                 "model with --model")
            faa = self._args.proteome
            fna = self._args.genome
            gff = self._args.gff
            model_path = self._args.model
            model_dir = os.path.dirname(model_path)

        ### SETTING UP DATA AND FUNCTIONS ###
        def calc_stoichiometry_df(df, counts):
            df["counts"] = pd.Series(counts)
            df["composition"] = df.counts / df.counts.sum()
            df["mass"] = df.mol_weight*df.composition
            df["stoichiometry"] = df.composition / df.mass.sum()
            df = df.fillna(0)
            return df

        DNA = pd.DataFrame({
                "id": ["C00131", "C00459", "C00286", "C00458"],
                "name": ["dATP", "dTTP", "dGTP", "dCTP"],
                "mol_weight": [0.491182, 0.482168, 0.507181, 0.467157]},
                index = ["A", "T", "G", "C"])

        RNA = pd.DataFrame({
                "id": ["C00002", "C00075", "C00044", "C00063"],
                "name": ["ATP", "UTP", "GTP", "CTP"],
                "mol_weight": [0.507181, 0.484141, 0.523180, 0.483156]},
                index = ["A", "U", "G", "C"])

        Prot = pd.DataFrame({
                "id":   ["C00041", "C00062", "C00152", "C00049",
                        "C00097", "C00025", "C00064", "C00037",
                        "C00135", "C00407", "C00123", "C00047",
                        "C00073", "C00079", "C00148", "C00065",
                        "C00188", "C00078", "C00082", "C00183"],
                "name": ["L-Alanine", "L-Arginine", "L-Asparagine",
                            "L-Aspartate", "L-Cysteine", "L-Glutamate",
                            "L-Glutamine", "L-Glycine", "L-Histidine",
                            "L-Isoleucine", "L-Leucine", "L-Lysine",
                            "L-Methionine", "L-Phenylalanine", "L-Proline",
                            "L-Serine", "L-Threonine", "L-Tryptophan",
                            "L-Tyrosine", "L-Valine"],
                "trna_ID":  ["C01635", "C01636", "C01637", "C01638",
                            "C01639", "C01641", "C01640", "C01642",
                            "C01643", "C01644", "C01645", "C01646",
                            "C01647", "C01648", "C01649", "C01650",
                            "C01651", "C01652", "C00787", "C01653"],
                "aatrna_ID":   ["C00886", "C02163", "C03402", "C02984",
                                "C03125", "C02987", "C02282", "C02412",
                                "C02988", "C03127", "C02047", "C01931",
                                "C02430", "C03511", "C02702", "C02553",
                                "C02992", "C03512", "C02839", "C02554"],
                "charge_rxn_ID":    ["R03038", "R03646", "R04212", "R05577",
                                    "R03650", "R05578", "R03905", "R03654",
                                    "R03655", "R03656", "R03657", "R03657",
                                    "R03659", "R03660", "R03661", "R03662",
                                    "R03663", "R03664", "R02918", "R03665"],
                #"charge_rxn_name":  ["L-Alanine:tRNA(Ala) ligase (AMP-forming)",
                #]
                "mol_weight":   [0.089093, 0.174201, 0.132118, 0.133103,
                                0.121158, 0.147129, 0.146145, 0.075067,
                                0.155155, 0.131173, 0.131173, 0.146188,
                                0.149211, 0.165189, 0.115131, 0.105093,
                                0.119119, 0.204225, 0.181189, 0.117146]},
                index = ["A", "R", "N", "D", "C", "E", "Q", "G",
                        "H", "I", "L", "K", "M", "F", "P", "S",
                        "T", "W", "Y", "V"])
        other = pd.DataFrame({
            "id":   ["C00020", "C00013", "C00035", "C00009",
                    "C00001", "C00008"],
            "name": ["AMP", "Diphosphate", "GDP", "Orthophosphate",
                    "H2O", "ADP"]
        })
        #template = pd.concat([Prot, DNA,
        #            RNA, other],
        #            ignore_index = True)
        #template.to_csv(sys.stdout)

        ### Handling the config file ###
        if self._args.generate_config:

            template = pd.concat([Prot[["id", "name"]], DNA[["id", "name"]],
                        RNA[["id", "name"]], other[["id", "name"]]],
                        ignore_index = True)
            template["custom_id"] = np.nan
            template.to_csv(sys.stdout, index = False)
            quit()

        if self._args.config:
            config = pd.read_csv(self._args.config)
            config.custom_id = config.custom_id.fillna(config.id)
            Prot = pd.merge(Prot, config, how = "inner", on = ["id", "name"])
            Prot.id = Prot.custom_id

        DNA_counts = Counter()
        RNA_counts = Counter()
        Prot_counts = Counter()

        yaml_args = {"default_flow_style": False,
                     "sort_keys": False,
                     "encoding": "utf-8",
                     "allow_unicode": True,
                     "width": float("inf")}

        ### Updating model.yaml ###
        logger.info("Updating model in {}".format(os.path.abspath(model_path)))

        with open(model_path, "r") as f:
            model_dict = yaml.safe_load(f)

        exists = False # Checking if biomass.yaml already included
        for file in model_dict["reactions"]:
            if "biomass.yaml" in file["include"]:
                exists = True
        if not exists:
            model_dict["reactions"].append({"include": "./biomass.yaml"})

        model_dict["biomass"] = self._args.biomass
        cell_compartment = model_dict["default_compartment"]
        if cell_compartment is None:
            cell_compartment = "c"
            model_dict["default_compartment"] = "c"

        os.rename(model_path, model_path + ".tmp")

        with open(os.path.join(model_path), "w") as f:
            yaml.dump(model_dict, f, **yaml_args)

        ### Check if needed compounds are in the model, download if not?
        needed_cpds = list(DNA.id) + list(RNA.id) + list(Prot.id) + \
                      list(Prot.aatrna_ID) + list(Prot.trna_ID) + \
                      list(other.id)
        model_reader = native.ModelReader.reader_from_path(model_path)
        compounds = [r.id for r in model_reader.parse_compounds()]
        missing_cpds = [cpd for cpd in needed_cpds if cpd not in compounds]
        if len(missing_cpds) > 0:
            logger.warning("The following compounds are required for biomass "
                            "equations but are missing from the model: "
                            "{}\n \033[1mEither create entries for these compounds "
                            "or specify corresponding names to existing "
                            "compounds using --config\033[0m".format(missing_cpds))

        ### DATA PROCESSING ###
        logger.info("Counting nucleotides and amino acids")

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
        ids_fwd = list(DNA.id) + ["C00002", "C00001"]
        stoich_fwd = list(DNA.stoichiometry) + [2*total, 2*total]
        compound_fwd = {Compound(id).in_compartment(cell_compartment):
                    round(-1*n, decimals) for id,n in zip(ids_fwd, stoich_fwd)}
        ids_rev = ["dna", "C00008", "C00009", "C00013"]
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
        ids_fwd = list(RNA.id) + ["C00001"]
        stoich_fwd = [RNA.loc["A"].stoichiometry + 2*total] +\
                    list(RNA.stoichiometry)[1:] + [2*total]
        compound_fwd = {Compound(id).in_compartment(cell_compartment):
                    round(-1*n, decimals) for id,n in zip(ids_fwd, stoich_fwd)}
        ids_rev = ["rna", "C00008", "C00009", "C00013"]
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
        ids_fwd = list(Prot.aatrna_ID) + ["C00044", "C00001"]
        stoich_fwd = list(Prot.stoichiometry) + [2*total, 2*total]
        compound_fwd = {Compound(id).in_compartment(cell_compartment):
                    round(-1*n, decimals) for id,n in zip(ids_fwd, stoich_fwd)}
        ids_rev = list(Prot.trna_ID) + ["protein", "C00035", "C00009"]
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
                        for id in ["dna", "rna", "protein"]}
        compound_rev = {Compound("Biomass").in_compartment(cell_compartment): 1}
        bio_rxn = Reaction(Direction.Forward, {**compound_fwd, **compound_rev})

        bio_rxn_entry = {"id": self._args.biomass,
                         "name": "Biomass",
                         "equation": str(bio_rxn),
                         "pathways": ["Biomass"]}

        ### Biomass sink reaction ###
        bio_sink = Reaction(Direction.Forward, {
                    Compound("Biomass").in_compartment(cell_compartment): -1})

        bio_sink_entry = {"id": "sink_biomass",
                         "name": "Biomass accumulation",
                         "equation": str(bio_sink),
                         "pathways": ["Biomass"]}

        ### tRNA charging reactions ###
        tRNA_rxns = []
        for index, row in Prot.iterrows():
            compound_fwd = {Compound(id).in_compartment(cell_compartment): -1
                            for id in ["C00002", row.id, row.trna_ID]}
            compound_rev = {Compound(id).in_compartment(cell_compartment): 1
                            for id in ["C00020", "C00013", row.aatrna_ID]}
            tRNA_rxn = Reaction(Direction.Both, {**compound_fwd,**compound_rev})
            tRNA_rxn_entry = {"id": row.charge_rxn_ID,
                              "name": "{} tRNA charging reaction".format(
                                                                row["name"]),
                              "equation": str(tRNA_rxn),
                              "pathways": ["Biomass"]}
            tRNA_rxns.append(tRNA_rxn_entry)

        ### dna, rna, protein, and Biomass compound formulas
        logger.info("Adding biomass compounds to {}/compounds.yaml".format(
                    os.path.abspath(model_dir)))
        cpd_entry = [{"id": id, "name": name, "charge": 0, "formula": "R"}
                        for id, name in [("protein", "Protein"),
                                         ("dna", "DNA"), ("rna", "RNA"),
                                         ("Biomass", "Biomass")]]
        with open(os.path.join(model_dir, "compounds.yaml"), "a") as f:
            yaml.dump(cpd_entry, f, **yaml_args)

        ### GENERATING biomass.yaml ###
        logger.info("Generating new biomass reactions in "
                    "{}/biomass.yaml".format(os.path.abspath(model_dir)))

        with open(os.path.join(model_dir, "biomass.yaml"), "w") as f:
            yaml.dump([dna_rxn_entry, rna_rxn_entry,
                       prot_rxn_entry, bio_rxn_entry,
                       bio_sink_entry] + tRNA_rxns,
                       f, **yaml_args)
