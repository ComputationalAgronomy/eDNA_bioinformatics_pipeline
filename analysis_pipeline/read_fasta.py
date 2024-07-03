from Bio import SeqIO

from analysis_pipeline.read_blast_csv import Reader

class ReadFasta(Reader):

    def __init__(self):
        super().__init__()
        self.seq_dict = {}

def read_fasta_file(seq_path, seq_type="Haplotype"):
    # type = "Haplotype"
    # type = "Amplicon"
    seq_dict = {}

    print(f"> {seq_type} FASTA files:  {seq_path}")
    n = 0

    fasta_sequences = SeqIO.parse(open(seq_path),'fasta')
    for fasta in fasta_sequences:
        name, seq = fasta.description, str(fasta.seq)
        if seq_type == "Amplicon":
            name = name.split(';')[0]
        seq_dict[name] = seq
        n += 1 # TODO(SW): if 'name' is unique, then this == len(seq_dict). If not unique, then this is the number of lines in the file, but do you need this number?

    # TODO(SW): You need to close() the file, or use with open() as file: ...
    print(f"{seq_type} Sequences:  {n} reads")

    return seq_dict