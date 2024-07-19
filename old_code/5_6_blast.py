import os

def blast_otu(in_dir, out_dir, db_path, otu_type='zotu', \
              specifiers='qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore', \
              cpu=2, maxhitnum=5, evalue=0.00001, percid=90, covperc=90):
    
    files = os.listdir(in_dir)
    number = []
    for filename in files:
        if filename.endswith(f'_{otu_type}.fasta'):
            number.append(filename.replace(f'_{otu_type}.fasta', ''))
    for num in number:
        cmd = f'blastn -db {db_path} -query {in_dir}{num}_{otu_type}.fasta -outfmt "10 {specifiers}" \
                -max_target_seqs {maxhitnum} -evalue {evalue} -perc_identity{percid} -qcov_hsp_perc {covperc} \
                -out {out_dir}{num}_{otu_type}.csv'
        cmd = cmd + f' -num_threads {cpu}' if cpu!=0 else cmd
        os.system(cmd)
