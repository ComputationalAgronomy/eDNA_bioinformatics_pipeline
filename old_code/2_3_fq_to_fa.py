import os 

def fq_to_fa(in_dir, out_dir, bbmap_dir):
    files = os.listdir(in_dir)
    number = []
    for filename in files:
        if filename.endswith('_cut.fastq'):
            number.append(filename.replace('_cut.fastq', ''))
    for num in number:
        cmd = f'bash {bbmap_dir}reformat.sh in={in_dir}{num}_cut.fastq out={out_dir}{num}_processed.fasta'
        os.system(cmd)
