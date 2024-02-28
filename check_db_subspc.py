from Bio import SeqIO
import os
import shutil
from pymsaviz import MsaViz

spc2seq = {}
subspc2seq = {}
species2sequence = {}
fasta_sequences = SeqIO.parse(open('./mifish_all_region.fasta'),'fasta')
for fasta in fasta_sequences:
    name, seq = fasta.description, str(fasta.seq)
    species = name.split('|')[-1]
    if len(species.split('_'))==1:
        species += '_sp.'

    if '_x_' in species:
        pass
    elif len(species.split('_'))==2:
        if species not in spc2seq:
            spc2seq[species]=[]
        spc2seq[species].append(seq)
    elif len(species.split('_'))>2:
        if species not in subspc2seq:
            subspc2seq[species]=[]
        subspc2seq[species].append(seq)
    else:
        print(name)

print(subspc2seq)

for species in spc2seq.keys():
    if species not in species2sequence:
        species2sequence[species] = {}
    species2sequence[species][species] = spc2seq[species]
    for subspecies in subspc2seq.keys():
        if species in subspecies:
            species2sequence[species][subspecies] = subspc2seq[subspecies]

# # put all sequence in one file to see alignment
# for species in species2sequence.keys():
#     num = len(species2sequence[species].keys())
#     if num>=2:
#         if not os.path.isdir(f'./check_db//{num}'):
#             os.makedirs(f'./check_db/{num}')
#         if not os.path.isdir('./tmp/'):
#             os.makedirs('./tmp/')
#         text = ''
#         file_string = ''
#         for spc, seq_list in species2sequence[species].items():    
#             umap_text = ''
#             for seq in seq_list:
#                 umap_text = umap_text + f'>{spc}\n{seq}\n'
#                 text = text + f'>{spc}\n{seq}\n'
#             with open(f'./tmp/{spc}.fasta', 'w')as file:
#                 file.write(umap_text)
#             file_string = file_string + f'./tmp/{spc}.fasta '

#         species = species.replace(':', '_')
#         species = species.replace('/', '_')
#         with open(f'./tmp/{species}_all.fasta', 'w')as file:
#             file.write(text)
#         cmd = f'usum {file_string} --maxdist 1.0 --termdist 1.0 --output ./check_db/{num}/{species} -f'
#         os.system(cmd)
#         cmd = f'clustalo -i ./tmp/{species}_all.fasta -o ./check_db/{num}/{species}/{species}.aln'
#         os.system(cmd)
#         mv = MsaViz(f'./check_db/{num}/{species}/{species}.aln', wrap_length=60, show_count=True, show_grid=False, show_consensus=True)
#         # mv.plotfig()
#         mv.savefig(f'./check_db/{num}/{species}/alignment.png')
#         shutil.rmtree('./tmp/')

# save every species, subspecies as independent file to plot umap
# for species in species2sequence.keys():
#     num = len(species2sequence[species].keys())
#     if num>=2:
#         if not os.path.isdir(f'./check_db/{num}'):
#             os.makedirs(f'./check_db/{num}')
#         if not os.path.isdir('./tmp/'):
#             os.makedirs('./tmp/')
#         file_string = ''
#         for spc, seq_list in species2sequence[species].items():
#             text = ''
#             for seq in seq_list:
#                 text = text + f'>{spc}\n{seq}\n'
#             spc = spc.replace(':', '_')
#             spc = spc.replace('/', '_')
#             with open(f'./tmp/{spc}.fasta', 'w') as file:
#                 file.write(text)
#             file_string += f'./tmp/{spc}.fasta '
#         cmd = f'usum {file_string}--maxdist 1.0 --termdist 1.0 --output ./check_db/{num}/{species} -f'
#         os.system(cmd)

