from Bio import SeqIO
import os


def check_sequence_number(in_dir, seq_fmt = 'fastq', dry = False):
    files = os.listdir(in_dir)
    total = 0
    for file in files:
        if file.endswith(seq_fmt) and dry==False:
            fastx_sequences = SeqIO.parse(open(check_dir+file), seq_fmt)
            for fastx in fastx_sequences:
                # title = fastx.description
                # size = int(title.split(';')[1].replace('size=',''))
                # total += size
                total += 1

    print(total)

def check_chimera_number(in_dir, prefix):
    files = os.listdir(in_dir)
    total = 0
    size_total = 0
    for file in files:
        if file.endswith(prefix):
            with open(check_dir+file, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if 'chfilter' in line and 'chimera' in line:
                        size = int(line.split(';')[1].replace('size=', ''))
                        size_total += size
                        total += 1
    print(total,"\n", size_total)

check_dir = './all_site/5_haploid/'
# check_sequence_number(in_dir=check_dir, seq_fmt='fasta', dry=False)
check_chimera_number(in_dir=check_dir, prefix='_zotu_report.txt')

# clean: 25632654
# merged: 11797855
# cut_primer: 11622521
# dereplicate: 515739(Unique sequences, size = 11622521)
# ZOTUs: # = 4246, size = 9873841)
# chimera: # = 1276, size = 43375