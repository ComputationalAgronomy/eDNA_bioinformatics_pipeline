from fasta import read_seq
from denoise_report import read_denoise_report
from blast_csv import read_hap_rank

class OneSample():

    def __init__(self, uniq_fasta_path, zotu_fasta_path, zotu_report_path, blast_csv_path):
        self.amp_seq = read_seq(uniq_fasta_path, seq_type="Amplicon")
        self.amp_size, self.hap2amp, self.hap_size = read_denoise_report(zotu_report_path)
        self.hap_seq = read_seq(zotu_fasta_path, seq_type="Haplotype")
        self.hap2rank = read_hap_rank(blast_csv_path)
