import subprocess
from os import listdir

path = './check_db/seq/10'
file_list = listdir(path)
file_list = [file.replace('.fasta', '') for file in file_list]
for filename in file_list:
    print(filename)
    cmd = f'blastn -db nt -remote -query {path}/{filename}.fasta -outfmt "10 qseqid sseqid stitle" -max_target_seqs 1 -out {path}/{filename}.csv'
    subprocess.run(cmd, shell=True)