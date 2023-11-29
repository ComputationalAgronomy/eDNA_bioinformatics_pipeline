from Bio import SeqIO
import copy

all = {}
complete = {}

fasta_sequences = SeqIO.parse(open('./mifish_all.fasta'),'fasta')
for fasta in fasta_sequences:
    name, seq = fasta.description, str(fasta.seq)
    all[name] = seq
print(len(all))
partial = copy.deepcopy(all)

fasta_sequences = SeqIO.parse(open('./mifish_complete.fasta'),'fasta')
for fasta in fasta_sequences:
    name, seq = fasta.description, str(fasta.seq)
    complete[name] = seq
print(len(complete))
rest = copy.deepcopy(complete)

for name, seq in all.items():
    repl = []
    for key, value in complete.items():
        if value == seq:
            repl.append(key)
    if len(repl)>1:
        print(repl)            
    

# print(rest)
# partial_list = [f'>{name}\n{seq}\n' for name, seq in partial.items()]
# text = ''.join(partial_list)

# with open('./mifish_partial.fasta', 'w') as file:
#     file.write(text)


# fasta_new = open('./database/mifish_complete_.fasta', 'w')
# fasta_sequences = SeqIO.parse(open('./database/mifish_complete.fasta'),'fasta')
# for fasta in fasta_sequences:
#     name, seq = fasta.description, str(fasta.seq)
#     name = name.replace(' ', '_')
#     fasta_new.write(f'>{name}\n{seq}\n')
# fasta_new.close()

# fasta_new = open('./mifish_complete_.fasta', 'w')
# fasta_sequences = SeqIO.parse(open('./mifish_complete.fasta'),'fasta')
# for fasta in fasta_sequences:
#     name, seq = fasta.description, str(fasta.seq)
#     name = name.replace(' ', '_')
#     fasta_new.write(f'>{name}\n{seq}\n')
# fasta_new.close()
