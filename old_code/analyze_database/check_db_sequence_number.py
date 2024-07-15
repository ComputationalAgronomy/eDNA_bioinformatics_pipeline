from Bio import SeqIO

total = 0
fasta_sequences = SeqIO.parse(open('./mifish_partial.fasta'),'fasta')
for fasta in fasta_sequences:
    total += 1

print(total)

# check complete genome
# species = []
# sub_species = []
# hybrid = []
# total = 0
# fasta_sequences = SeqIO.parse(open('./mifish_complete.fasta'),'fasta')
# for fasta in fasta_sequences:
#     total += 1
#     name, seq = fasta.description, str(fasta.seq)
#     spc = name.split('|')[-1]
#     spc_list = spc.split('_')
#     if len(spc_list) == 1:
#         spc += '_sp.'
    
#     if '_x_' in spc and spc not in hybrid:
#         hybrid.append(spc)
#         continue

#     if 'sp.' in spc and spc not in species:
#         species.append(spc)
#         continue

#     if len(spc_list) == 2 and spc not in species:
#         species.append(spc)
#         continue
    
#     if len(spc_list) == 3 and spc not in sub_species:
#         sub_species.append(spc)
#         continue

# print(total)
# print(len(species))
# print(len(sub_species))
# print(len(hybrid))