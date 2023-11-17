import os


def cluster_otu(in_dir, out_dir, minsize=2, cpu=2):
    files = os.listdir(in_dir)
    number = []
    for filename in files:
        if filename.endswith('_derep.fasta'):
            number.append(filename.replace('_derep.fasta', ''))
    for num in number:
        cmd = f'usearch -cluster_otus {in_dir}{num}_derep.fasta -minsize {minsize} -threads {cpu} \
                -otus {out_dir}{num}_otu.fasta -uparseout {out_dir}{num}_otu_report.txt -relabel Otu \
                >{out_dir}{num}_report.txt 2>&1' 
        os.system(cmd)
