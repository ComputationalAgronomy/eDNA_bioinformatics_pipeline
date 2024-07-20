import os

def dereplicate(in_dir, out_dir, cpu=2):
    files = os.listdir(in_dir)
    number = []
    for filename in files:
        if filename.endswith('_processed.fasta'):
            number.append(filename.replace('_processed.fasta', ''))
    for num in number:
        cmd = f'usearch -fastx_uniques {in_dir}{num}_processed.fasta -threads {cpu} \
                -sizeout -relabel Uniq -fastaout {out_dir}{num}_derep.fasta \
                >{out_dir}{num}_report.txt 2>&1' 
        os.system(cmd)
