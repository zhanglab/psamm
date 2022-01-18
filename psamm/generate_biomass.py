from collections import Counter
from six import add_metaclass, iteritems, itervalues, text_type, PY3
import abc
from psamm.datasource.reaction import Reaction, Compound, Direction
try:
    from Bio import SeqIO
except ModuleNotFoundError:
    quit("No biopython package found. Please run <pip install biopython>")


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
            logger.debug('Command failure caused by exception!', exc_info=exc)
        sys.exit(1)


class InputError(Exception):
    """Exception used to signal a general input error."""


class main(Command):
    """Generate a database of compounds and reactions"""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument('--genome', metavar='path',
            help='Path to the genome in fasta format')
        parser.add_argument('--proteome', metavar='path',
            help='Path to the proteome in fasta format')
        parser.add_argument('--gff', metavar = 'path',
            help='Specify path to gff containing transcript annotations. '
                 'Only annotations specified as \'CDS\' in the third column will be used. '
                 'Annotations must correspond to sequences in the --genome file.')
        parser.add_argument('--model', metavar='path',
            help='Path to the model directory')
        ## option for name of biomass function
        ## option for renaming biomass components based on table?
        ## update model.yaml
        ## name your new branch genBiomass
        super(main, cls).init_parser(parser)

    def run(self):
        print("HELLOO")
        """Entry point for the biomass reaction generation script"""
        if not self._args.genome:
            raise InputError('Please specify a path to the genome in fasta format')
        if not self._args.proteome:
            raise InputError('Specify path to proteome in fasta format')
        if not self._args.gff:
            raise InputError('Specify path to gff containing transcript annotations. '
                 'Only annotations specified as \'CDS\' in the third column will be used. '
                 'Annotations must correspond to sequences in the --genome file.')
        if not self._args.model:
            raise InputError('Please specify a directory to an existing model')
        faa = self._args.proteome
        fna = self._args.genome
        gff = self._args.gff

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

        # Parsing through the gff and counting bases in coding regions to approximate
        # the total RNA base counts
        with open(gff, "r") as gff_file:
            for line in gff_file:
                if not line.startswith("#"):
                    l = line.strip().split("\t")
                    if l[2] == "CDS":
                        try:
                            seq_entry = genome[l[0]]
                        except KeyError:
                            print("WARNING: an annotation in gff file does not match any\
                            sequences in fasta file, skipping...")
                            continue
                        ## FIX STRANDEDNESS:: if + use main strand, if - use reverse complement
                        RNA_counts += Counter(seq_entry.seq[int(l[3]):int(l[4])])

        RNA_counts["U"] = RNA_counts.pop("T") #renames T to U for RNA
        def calc_stoichiometry(counts, weights):
            total_count = sum(counts.values())
            composition = {base: counts[base]/total_count for base in counts.keys()}
            mass = sum([composition[base]*weights[base] for base in counts.keys()])
            stoichiometry = {base: round(composition[base]/mass, 6) for base in counts.keys()}
            return stoichiometry

        DNA_stoichiometry = calc_stoichiometry(DNA_counts, {"A": 0.4911816, "T": 0.4821683, "G": 0.507181, "C": 0.4671569})
        RNA_stoichiometry = calc_stoichiometry(RNA_counts, {"A": 0.507181, "U": 0.484141, "G": 0.52318, "C": 0.483156})
        Prot_stoichiometry = calc_stoichiometry(Prot_counts, {"A": 0.089093, "R": 0.174201, "N": 0.132118, "D": 0.133103, "C": 0.121158, "E": 0.147129, "Q": 0.146145, "G": 0.075067, "H": 0.155155, "I": 0.131173, "L": 0.131173, "K": 0.146188, "M": 0.149211, "F": 0.165189, "P": 0.115131, "S": 0.105093, "T": 0.119119, "W": 0.204225, "Y": 0.181189, "V": 0.117146})
        print(DNA_stoichiometry, RNA_stoichiometry, Prot_stoichiometry)




        #total_bases = sum(DNA_counts.values())
        #DNA_base_composition = {base: DNA_counts[base]/total_bases for base in DNA_counts.keys()}
        #DNA_weights = {"A": 0.4911816, "T": 0.4821683, "G": 0.507181, "C": 0.4671569}
        #DNA_mass = sum([DNA_base_composition[base]*DNA_weights[base] for base in DNA_counts.keys()])
        #DNA_stoichiometry = {base: round(DNA_base_composition[base]/DNA_mass, 6) for base in DNA_counts.keys()}
#
        #rna_bases = sum(RNA_counts.values())
        #RNA_base_composition = {base: RNA_counts[base]/rna_bases for base in RNA_counts.keys()}
        #RNA_weights = {"A": 0.507181, "U": 0.484141, "G": 0.52318, "C": 0.483156}
        #RNA_mass = sum([RNA_base_composition[base]*RNA_weights[base] for base in RNA_counts.keys()])
        #RNA_stoichiometry = {base: round(RNA_base_composition[base]/RNA_mass, 6) for base in RNA_counts.keys()}
#
        #total_amino_acids = sum(Prot_counts.values())
        #Prot_base_composition = {aa: Prot_counts[aa]/total_amino_acids for aa in Prot_counts.keys()}
        #Prot_weights = {"A": 0.089093, "R": 0.174201, "N": 0.132118, "D": 0.133103, "C": 0.121158, "E": 0.147129, "Q": 0.146145, "G": 0.075067, "H": 0.155155, "I": 0.131173, "L": 0.131173, "K": 0.146188, "M": 0.149211, "F": 0.165189, "P": 0.115131, "S": 0.105093, "T": 0.119119, "W": 0.204225, "Y": 0.181189, "V": 0.117146}
        #Prot_mass = sum([Prot_base_composition[aa]*Prot_weights[aa] for aa in Prot_counts.keys()])
        #Prot_stoichiometry = {aa: round(Prot_base_composition[aa]/Prot_mass, 6) for aa in Prot_counts.keys()}
        #print(DNA_stoichiometry, RNA_stoichiometry, Prot_stoichiometry)
