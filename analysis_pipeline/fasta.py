from Bio import SeqIO

def read_seq(seq_path, seq_type="Haplotype"):
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
        n += 1
    
    print(f"{seq_type} Sequences:  {n} reads")

    return seq_dict