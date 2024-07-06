from Bio import SeqIO

from read_blast_csv import Reader

class FastaReader(Reader):

    def __init__(self):
        super().__init__()
        self.seq_dict = {}

    @staticmethod
    def parse_fasta(seq_path: str, seq_type: str):
        """
        Parse the fasta file and return a list of tuples containing sequence names and sequences.

        :param seq_path: The path to the fasta file.
        :param seq_type: The type of sequence, either "Haplotype" or "Amplicon".
        :return: A list of tuples (name, seq).
        """
        sequences = []
        with open(seq_path) as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                name, seq = record.description, str(record.seq)
                if seq_type == "Amplicon":
                    name = name.split(';')[0]
                sequences.append((name, seq))

        return sequences

    def update_seq_dict(self, sequences):
        """
        Update the dictionary 'self.seq_dict' with the provided sequences.

        :param sequences: A list of tuples (name, seq).
        """
        for name, seq in sequences:
            self.seq_dict[name] = seq

    def read_fasta(self, seq_path: str, seq_type: str ="Haplotype"):
        """
        Read a fasta file and update a dictionary of sequences.

        :param seq_path: The path to the fasta file.
        :param seq_type: The type of sequence,  either "Haplotype" or "Amplicon", default is "Haplotype".
        """
        print(f"> {seq_type} FASTA files:  {seq_path}")

        sequences = FastaReader.parse_fasta(seq_path, seq_type)
        self.update_seq_dict(sequences)

        read_count = len(self.seq_dict)
        print(f"{seq_type} Sequences:  {read_count} reads")