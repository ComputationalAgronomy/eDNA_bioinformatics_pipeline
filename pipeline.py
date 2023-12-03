import os
import pandas as pd

def get_prefix_with_suffix(in_dir, suffix):
    files = os.listdir(in_dir)
    prefix = []
    for filename in files:
        if filename.endswith(suffix):
            prefix.append(filename.replace(suffix, ''))
    return prefix

def merge_fq(in_dir, out_dir, maxdiff=5, pctid=90, cpu=2):
    prefix = get_prefix_with_suffix(in_dir, '_R1.fastq')
    for num in prefix:
        cmd = f'usearch -fastq_mergepairs {in_dir}{num}_R1.fastq -fastqout {out_dir}{num}_merge.fastq -threads {cpu} \
                -fastq_maxdiffs {maxdiff} -fastq_pctid {pctid} -report {out_dir}{num}_report.txt'
        os.system(cmd)

def cut_adapt(in_dir, out_dir, rm_p_5='GTCGGTAAAACTCGTGCCAGC', rm_p_3='CAAACTGGGATTAGATACCCCACTATG', min_read_len=204, max_read_len=254, cpu=2):
    prefix = get_prefix_with_suffix(in_dir, '_merge.fastq')
    for num in prefix:
        os.system(f'cutadapt -g "{rm_p_5};max_error_rate=0.15...{rm_p_3};max_error_rate=0.15" -j {cpu} \
                    {in_dir}{num}_merge.fastq --discard-untrimmed \
                    -m {min_read_len-len(rm_p_5)-len(rm_p_3)} -M {max_read_len-len(rm_p_5)-len(rm_p_3)} \
                    >{out_dir}{num}_cut.fastq 2>{out_dir}{num}_report.txt')

def fq_to_fa(in_dir, out_dir, bbmap_dir):
    prefix = get_prefix_with_suffix(in_dir, '_cut.fastq')
    for num in prefix:
        cmd = f'bash {bbmap_dir}reformat.sh in={in_dir}{num}_cut.fastq out={out_dir}{num}_processed.fasta'
        os.system(cmd)

def dereplicate(in_dir, out_dir, cpu=2):
    prefix = get_prefix_with_suffix(in_dir, '_processed.fasta')
    for num in prefix:
        cmd = f'usearch -fastx_uniques {in_dir}{num}_processed.fasta -threads {cpu} \
                -sizeout -relabel Uniq -fastaout {out_dir}{num}_derep.fasta \
                >{out_dir}{num}_report.txt 2>&1'
        os.system(cmd)

def cluster_otu(in_dir, out_dir, minsize=2, cpu=2):
    prefix = get_prefix_with_suffix(in_dir, '_derep.fasta')
    for num in prefix:
        cmd = f'usearch -cluster_otus {in_dir}{num}_derep.fasta -minsize {minsize} -threads {cpu} \
                -otus {out_dir}{num}_otu.fasta -uparseout {out_dir}{num}_otu_report.txt -relabel Otu \
                >{out_dir}{num}_report.txt 2>&1'
        os.system(cmd)

def cluster_zotu(in_dir, out_dir, minsize=8, cpu=2):
    prefix = get_prefix_with_suffix(in_dir, '_derep.fasta')
    for num in prefix:
        cmd = f'usearch -unoise3 {in_dir}{num}_derep.fasta -minsize {minsize} -threads {cpu} \
                -zotus {out_dir}{num}_zotu.fasta -tabbedout {out_dir}{num}_zotu_report.txt \
                >{out_dir}{num}_report.txt 2>&1'
        os.system(cmd)

def blast_otu(in_dir, out_dir, db_path, otu_type='otu', cpu=0, maxhitnum=5, specifiers='qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'):
    if otu_type in ['otu', 'zotu']:
        pass
    else:
        print("Error: 'otu_type' must be 'otu' or 'zotu'.")
        return

    prefix = get_prefix_with_suffix(in_dir, f'_{otu_type}.fasta')
    for num in prefix:
        cmd = f'blastn -db {db_path} -query {in_dir}{num}_{otu_type}.fasta -outfmt "10 {specifiers}" \
                -max_target_seqs {maxhitnum} -evalue 0.00001 -qcov_hsp_perc 90 \
                -out {out_dir}{num}_{otu_type}.csv'
        cmd = cmd + f' -num_threads {cpu}' if cpu!=0 else cmd
        os.system(cmd)

    for num in prefix:
        otu_hit = pd.read_csv(f'{out_dir}{num}_{otu_type}.csv', header=None)
        species_name = [x[2] for x in otu_hit[1].str.split('|')]
        otu_hit[1] = [x[1] for x in otu_hit[1].str.split('|')]
        otu_hit.insert(2, '2', species_name)
        otu_hit.to_csv(f'{out_dir}{num}_{otu_type}.csv', index=False, header=None)

if __name__ == '__main__':

    in_dir = './cleandata/5_haploid/otu/'
    out_dir = './cleandata/6_blastn/mifish_db/'

    # merge_fq(in_dir=in_dir, out_dir=out_dir, cpu=6)

    # cut_adapt(in_dir=in_dir, out_dir=out_dir, cpu=6)

    # bbmap_dir = './bbmap/'
    # fq_to_fa(bbmap_dir=bbmap_dir, in_dir=in_dir, out_dir=out_dir)

    # dereplicate(in_dir=in_dir, out_dir=out_dir, cpu=6)

    # cluster_otu(in_dir=in_dir, out_dir=out_dir, minsize=8, cpu=16)
    # cluster_zotu(in_dir=in_dir, out_dir=out_dir, minsize=8, cpu=16)

    # mifish_path = './database/mifishdb'
    # blast_otu(in_dir=in_dir, out_dir=out_dir, db_path=mifish_path, maxhitnum=1, otu_type='otu', cpu=16)

    # # ncbi_path = 'nt -remote'
    # blast_otu(in_dir=in_dir, out_dir=out_dir, db_path=ncbi_path, out_name='ncbi', otu_type='zotu')
    # read_blast_csv(in_dir=in_dir, out_dir=out_dir)