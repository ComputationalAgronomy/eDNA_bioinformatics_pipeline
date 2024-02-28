from Bio import SeqIO
import re

match_pattern = re.compile(r"GTCGGTAAAACTCGTGCCAGC[ATCG]+CAAACTGGGATTAGATACCCCACTATG")

fasta_new = open('./mifish_all_region.fasta', 'w')
fasta_sequences = SeqIO.parse(open('./mifish_all.fasta'),'fasta')
for fasta in fasta_sequences:
    name, seq = fasta.description, str(fasta.seq)
    match = match_pattern.search(seq)
    if match:
        seq_mifish = match.group()
        seq_mifish = seq_mifish[21:-27]
        fasta_new.write(f'>{name}\n{seq_mifish}\n')
fasta_new.close()

