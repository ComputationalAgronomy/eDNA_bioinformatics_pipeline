import csv
from Bio import SeqIO
import subprocess
import os
import shutil
from pymsaviz import MsaViz

family2genus = {}
with open('./database/lineage.csv') as in_handle:
    csvReader = csv.DictReader(in_handle)
    for row in csvReader:
        if row['family_name'] not in family2genus:
            family2genus[row['family_name']] = []
        family2genus[row['family_name']].append(row['genus_name'])

genus2species = {}
species2seq = {}
fasta_sequences = SeqIO.parse(open('./mifish_12s.fasta'),'fasta')
for fasta in fasta_sequences:
    species = fasta.description.split('|')[-1]
    genus = species.split('_')[0]
    if genus not in genus2species:
        genus2species[genus] = []
    if species not in genus2species[genus]:
        genus2species[genus].append(species)
    seq = fasta.seq
    if species not in species2seq:
        species2seq[species] = []
    species2seq[species].append(seq)


# for genus, species_list in genus2species.items():
#     file_string = ''
#     align_seq = ''
#     if not os.path.isdir('./tmp/'):
#         os.makedirs('./tmp/')
#     for species in species_list:
#         umap_seq = ''
#         for seq in species2seq[species]:
#             umap_seq += f'>{species}\n{seq}\n'
#             align_seq += f'>{species}\n{seq}\n'
#         with open(f'./tmp/{species}.fa', 'w')as file:
#             file.write(umap_seq)
#         file_string += f'./tmp/{species}.fa '
#     with open(f'./tmp/seq.fa', 'w')as file:
#         file.write(align_seq)

#     cmd = f'usum {file_string} --maxdist 1.0 --termdist 1.0 --output ./check_db/genus/{genus} -f'
#     subprocess.run(cmd, shell=True)
#     cmd = f'clustalo -i ./tmp/seq.fa -o ./check_db/genus/{genus}/seq.aln --distmat-out=./check_db/genus/{genus}/seq_dist.mat --guidetree-out=./check_db/genus/{genus}/seq_tree.dnd --full --force'
#     subprocess.run(cmd, shell=True)
#     if os.path.exists(f'./check_db/genus/{genus}/seq.aln'):
#         mv = MsaViz(f'./check_db/genus/{genus}/seq.aln', wrap_length=60, show_count=True, show_grid=False, show_consensus=True)
#         mv.savefig(f'./check_db/genus/{genus}/alignment.png')

#     shutil.rmtree('./tmp/')

for family, genus_list in family2genus.items():
    file_string = ''
    align_seq = ''
    if not os.path.isdir('./tmp2/'):
        os.makedirs('./tmp2/')
    for genus in genus_list:
        if genus in genus2species:
            umap_seq = ''
            species_list = genus2species[genus]
            for species in species_list:
                align_seq += f'>{species}\n{species2seq[species]}\n'
                umap_seq += f'>{species}\n{species2seq[species]}\n'
            with open(f'./tmp2/{genus}.fa', 'w')as file:
                file.write(umap_seq)
            file_string += f'./tmp2/{genus}.fa '         
    with open(f'./tmp2/seq.fa', 'w')as file:
        file.write(align_seq)

    cmd = f'usum {file_string} --maxdist 1.0 --termdist 1.0 --output ./check_db/family/{family} -f'
    subprocess.run(cmd, shell=True)

    cmd = f'clustalo -i ./tmp2/seq.fa -o ./check_db/family/{family}/seq.aln --distmat-out=./check_db/family/{family}/seq_dist.mat --guidetree-out=./check_db/family/{family}/seq_tree.dnd --full --force'
    subprocess.run(cmd, shell=True)
    if os.path.exists(f'./check_db/family/{family}/seq.aln'):
        mv = MsaViz(f'./check_db/family/{family}/seq.aln', wrap_length=60, show_count=True, show_grid=False, show_consensus=True)
        mv.savefig(f'./check_db/family/{family}/alignment.png')

    shutil.rmtree('./tmp2/')