from Bio import SeqIO

complete_partial = {}
complete = {}
fasta_sequences = SeqIO.parse(open('./database/mifish_complete_partial.fasta'),'fasta')
for fasta in fasta_sequences:
    name, seq = fasta.description, str(fasta.seq)
    complete_partial[name] = seq
 
fasta_sequences = SeqIO.parse(open('./database/mifish_complete.fasta'),'fasta')
for fasta in fasta_sequences:
    name, seq = fasta.description, str(fasta.seq)
    complete[name] = seq

partial = set(complete_partial) - set(complete)

print("found difference.")

partial_list = [f'>{k}\n{complete_partial[k]}' for k in partial]


text = f'\n'.join(partial_list)
print('made text.')
# text = ""
# for k in partial:
#     print(k)
#     # text += f'>{k}\n{complete_partial[k]}\n'

print("start writing.")
with open('./database/mifish_partial.fasta', 'w') as file:
    file.write(text)