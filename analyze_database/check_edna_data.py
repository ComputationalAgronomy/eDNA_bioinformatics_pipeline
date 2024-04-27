from Bio import SeqIO
import os
import pandas as pd
import re


def check_sequence_number(in_dir, prefix = '', seq_fmt = 'fastq', dry = False):
    files = os.listdir(in_dir)
    total = 0
    total_size= 0
    for file in files:
        if file.endswith(prefix+seq_fmt) and dry==False:
            fastx_sequences = SeqIO.parse(open(in_dir+file), seq_fmt)
            for fastx in fastx_sequences:
                title = fastx.description
                size = int(title.split(';')[1].replace('size=',''))
                total_size += size
                total += 1

    print(total)
    print(total_size)

def check_denoise(in_dir, prefix):
    zotu = 0
    zotu_size = 0
    chimera = 0
    chimera_size = 0
    zotu8 = 0
    zotu8_size = 0
    chimera8 = 0
    chimera8_size = 0

    files = os.listdir(in_dir)
    for file in files:
        if file.endswith(prefix):
            with open(in_dir+file, 'r') as f:
                amp_size = {}
                zotu_num = 1
                chi_num = 1
                for line in f.readlines():
                    line_list = line.split(';')
                    haploid= line_list[0]
                    if 'denoise' in line:
                        size = int(line_list[1].replace('size=','')) 
                        if 'amp' in line:
                            amp = haploid
                        elif 'shifted' in line_list[2] or 'bad' in line_list[2]:
                            amp = re.search(r'top=Uniq\d+', line_list[3]).group(0).replace('top=','')

                        if amp not in amp_size:
                            amp_size[amp] = 0
                        amp_size[amp] += size

                    elif 'chfilter' in line_list[2]:
                        if 'zotu' in line_list[2]:
                            k_new = f'Zotu{zotu_num}'
                            amp_size[k_new] = amp_size.pop(haploid)
                            zotu_num += 1
                        elif 'chimera' in line_list[2]:
                            k_new = f'Chimera{chi_num}'
                            amp_size[k_new] = amp_size.pop(haploid)
                            chi_num += 1

    
                for amp, size in amp_size.items():
                    if 'Zotu' in amp:
                        zotu += 1
                        zotu_size += size
                        if size >= 8:
                            zotu8 += 1
                            zotu8_size += size
                    elif 'Chimera' in amp:
                        chimera += 1
                        chimera_size += amp_size[amp]
                        if size >= 8:
                            chimera8 += 1
                            chimera8_size += size
    print('zotu:', zotu)
    print('zotu_size:', zotu_size)
    print('chimera:', chimera)
    print('chimera_size:', chimera_size)
    print('zotu8:', zotu8)
    print('zotu8_size:', zotu8_size)
    print('chimera8:', chimera8)
    print('chimera8_size:', chimera8_size)

def check_hit_number(in_dir, prefix):
    files = os.listdir(in_dir)
    total = 0

    for file in files:
        if file.endswith(prefix):
            results = pd.read_csv(in_dir+file)
            total += len(results)
    
    print(total)

check_dir = './data/all_site/5_haploid/'
# check_sequence_number(in_dir=check_dir, prefix='2303-H03_derep.', seq_fmt='fasta', dry=False)
check_denoise(in_dir=check_dir, prefix='_zotu_report.txt')
# check_hit_number(in_dir=check_dir, prefix = '.csv')

# clean: 25632654
# merged: 11797855
# cut_primer: 11622521
# dereplicate: 515739(Unique sequences, size = 11622521)
# ZOTUs: = 4246, size = 10895700
# chimera: = 1276, size = 44485
# low abundance ZOTU removed: 20202 (640902)
# low abundance chimera removed: 21450 (41434)