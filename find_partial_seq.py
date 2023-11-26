from Bio import SeqIO

complete_partial = {}
complete = {}
fasta_sequences = SeqIO.parse(open('./database/mifish_complete_partial.fasta'),'fasta')
for fasta in fasta_sequences:
    name, seq = fasta.description, str(fasta.seq)
    complete_partial[name] = seq
print(len(complete_partial))
partial = complete_partial.copy()

fasta_sequences = SeqIO.parse(open('./database/mifish_complete.fasta'),'fasta')
for fasta in fasta_sequences:
    name, seq = fasta.description, str(fasta.seq)
    complete[name] = seq
print(len(complete))

for name, seq in complete_partial.items():
    if seq in complete.values():
        partial.pop(name)
print(len(partial))

partial_list = [f'>{name}\n{seq}' for name, seq in partial.items()]
text = ''.join(f"{row}\n" for row in partial_list)

with open('./database/mifish_partial.fasta', 'w') as file:
    file.write(text)
