import os
from Bio import SeqIO
import re

genome_list = [file for file in os.listdir('./mitogenomes')]
annotation_list = [file for file in os.listdir('./mitoannotations') if file.endswith('.txt') and not file.endswith('.NCBI.txt')]
for i in range(len(genome_list)):
    with open(f'./mitoannotations/{annotation_list[i]}', 'r')as file:
        lines = file.readlines()
        for line in lines:
            if '12S rRNA' in line:
                region = re.findall(r'\d+..\d+', line)[0]
                continue
        print(region)
        start, end = region.split('..')

    text = ''
    for fasta in SeqIO.parse(f'./mitogenomes/{genome_list[i]}','fasta'):
        species_name = fasta.description.split('|')[-1]
        seq = fasta.seq[int(start):int(end)+1]
        text += f'>{fasta.description}\n{seq}\n'
    with open(f'./mitoannotations_12s/{genome_list[i]}', 'w')as file:
        file.write(text)
