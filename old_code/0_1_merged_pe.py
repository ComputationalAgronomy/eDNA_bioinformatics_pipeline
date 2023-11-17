import os

def merge_fastq(fastq_dir, save_dir, cpu=2, maxdiff=5, pctid=90):
    files = os.listdir(fastq_dir)
    number = []
    for filename in files:
        if filename.endswith('_R1.fastq'):
            number.append(filename.replace('_R1.fastq', ''))
    for num in number:
        cmd = f'usearch -fastq_mergepairs {fastq_dir}{num}_R1.fastq -fastqout {save_dir}{num}_merged.fastq -threads {cpu} \
                -fastq_maxdiffs {maxdiff} -fastq_pctid {pctid} -report {save_dir}{num}_report.txt' 
        os.system(cmd)
