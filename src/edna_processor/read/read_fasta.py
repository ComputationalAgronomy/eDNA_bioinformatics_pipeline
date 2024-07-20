from Bio import SeqIO

from edna_processor.read.read_blast_csv import Reader
from edna_processor.utils.base_logger import logger

class FastaReader(Reader):

    def __init__(self):
        super().__init__()
        self.seq_dict = {}

    @staticmethod
    def parse_fasta(seq_path: str, seq_type: str) -> dict[str, str]:
        """
        Parse the fasta file and return the dictionary 'seq_dict' containing sequence names and sequences.

        :param seq_path: The path to the fasta file.
        :param seq_type: The type of sequence, either "Haplotype" or "Amplicon".
        :return: A dictionary containing sequence names and sequences.
        """
        seq_dict = {}
        with open(seq_path) as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                name, seq = record.description, str(record.seq)
                if seq_type == "Amplicon":
                    name = name.split(';')[0]
                seq_dict[name] = seq
        return seq_dict

    def read_fasta(self, seq_path: str, seq_type: str ="Haplotype") -> None:
        """
        Read a fasta file and update the dictionary 'seq_dict' with sequence names and sequences.

        :param seq_path: The path to the fasta file.
        :param seq_type: The type of sequence,  either "Haplotype" or "Amplicon". Default is "Haplotype".
        """
        logger.info(f"Reading {seq_type} FASTA files:  {seq_path}")

        self.seq_dict = self.parse_fasta(seq_path, seq_type)

        read_count = len(self.seq_dict)
        logger.info(f"Read finished. {seq_type} Sequences:  {read_count} reads")