from Bio import SeqIO

def split_lines(seq):
        lines = [seq[i:i+59] for i in range(0, len(seq), 60)]
        return '\n'.join(lines)

text = ''
fasta_sequences = SeqIO.parse(open('./mifish_all.fasta'),'fasta')
for fasta in fasta_sequences:
    name, seq = fasta.description.split('|')[-1], str(fasta.seq)
    if 'Mugil_' in name and seq.startswith('CACCGCGG'):
        text += f'>{fasta.description}\n{seq}\n'

print(text)