import os


def cluster_zotu(in_dir, out_dir, minsize=8, cpu=2):
    files = os.listdir(in_dir)
    number = []
    for filename in files:
        if filename.endswith('_derep.fasta'):
            number.append(filename.replace('_derep.fasta', ''))
    for num in number:
        cmd = f'usearch -unoise3 {in_dir}{num}_derep.fasta -minsize {minsize} -threads {cpu} \
                -zotus {out_dir}{num}_zotu.fasta -tabbedout {out_dir}{num}_zotu_report.txt \
                >{out_dir}{num}_report.txt 2>&1'
        os.system(cmd)
