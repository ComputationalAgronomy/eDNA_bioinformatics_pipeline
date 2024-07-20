from Bio import SeqIO
import subprocess
from pymsaviz import MsaViz

complete = {}
fasta_sequences = SeqIO.parse(open('./mifish_complete.fasta'),'fasta')
for fasta in fasta_sequences:
    name, seq = fasta.description, str(fasta.seq)
    complete[name] = seq

rev_dict = {}
for key, value in complete.items():
    rev_dict.setdefault(value, set()).add(key)

dupl_num = 1
for values in rev_dict.values():
    if len(values)>1:
        fasta_new = open(f'./dupl_{dupl_num}.fasta', 'w')
        for name in values:
            fasta_new.write(f'>{name}\n{complete[name]}\n')
        fasta_new.close()

        cmd = f'clustalo -i ./dupl_{dupl_num}.fasta -o ./dupl_{dupl_num}.aln'
        subprocess.run(cmd, shell=True)
        mv = MsaViz(f'./dupl_{dupl_num}.aln', wrap_length=60, show_count=True, show_grid=False, show_consensus=True)
        mv.savefig(f'./dupl_{dupl_num}_align.png')
        dupl_num += 1
