from fasta import read_seq
from denoise_report import read_denoise_report
from blast_csv import read_hap_rank
import os
import time

class OneSample():

    def __init__(self, uniq_fasta_path, zotu_fasta_path, zotu_report_path, blast_csv_path):
        self.amp_seq = read_seq(uniq_fasta_path, seq_type="Amplicon")
        self.amp_size, self.hap2amp, self.hap_size = read_denoise_report(zotu_report_path)
        self.hap_seq = read_seq(zotu_fasta_path, seq_type="Haplotype")
        self.hap2rank = read_hap_rank(blast_csv_path)

class IntegrateSamples():

    def __init__(self, load_path=None):
        self.sample_data = {}
        self.sample_id_list = []

        if load_path:
            self.load_sample_data(load_path)

    def import_samples(self, read_dir, sample_id_list=None): # read_dir: path to the directory containing the "4_derep", "5_denoise", "6_blast" folders.

        start_time = time.time()

        print(f"> Reading samples from: {read_dir}.")
        if os.path.exists(os.path.join(read_dir, "4_derep"))==False:
            print("> Error: 4_derep folder does not exist in the specified directory.")
            return
        if os.path.exists(os.path.join(read_dir, "5_denoise"))==False:
            print("> Error: 5_denoise folder does not exist in the specified directory.")
            return
        if os.path.exists(os.path.join(read_dir, "6_blast"))==False:
            print("> Error: 6_blast folder does not exist in the specified directory.")
            return

        if sample_id_list is None:
            print("> No sample id list provided, reading sample id list from the directory.")
            file_list = os.listdir(f'{read_dir}/4_derep')
            sample_id_list = [file.replace('_uniq.fasta', '') for file in file_list if file.endswith('_uniq.fasta')]
        else:
            print("> Specified sample id list. Start checking whether files exist.")
            for sample_id in sample_id_list:
                if os.path.exists(os.path.join(read_dir, "4_derep", f"{sample_id}_uniq.fasta"))==False:
                    print(f"> Error: {sample_id}_uniq.fasta does not exist in the specified directory.")
                    return
                if os.path.exists(os.path.join(read_dir, "5_denoise", f"{sample_id}_zotu.fasta"))==False:
                    print(f"> Error: {sample_id}_zotu.fasta does not exist in the specified directory.")
                    return
                if os.path.exists(os.path.join(read_dir, "6_blast", f"{sample_id}_blast.csv"))==False:
                    print(f"> Error: {sample_id}_blast.csv does not exist in the specified directory.")
                    return
            print("> All files exist.")

            check_file(read_dir, sample_id_list)

        self.sample_id_list.extend(sample_id_list)

        for sample_id in sample_id_list:

            uniq_fasta_path = os.path.join(read_dir, "4_derep", f"{sample_id}_uniq.fasta")
            zotu_fasta_path = os.path.join(read_dir, "5_denoise", f"{sample_id}_zotu.fasta")
            zotu_report_path = os.path.join(read_dir, "5_denoise", f"{sample_id}_zotu_report.txt")
            blast_csv_path = os.path.join(read_dir, "6_blast", f"{sample_id}_blast.csv")

            self.sample_data[sample_id] = OneSample(uniq_fasta_path=uniq_fasta_path, zotu_fasta_path=zotu_fasta_path, zotu_report_path=zotu_report_path, blast_csv_path=blast_csv_path)

        seconds = time.time() - start_time
        print(f"> {len(sample_id_list)} Samples read.")
        print("Time Taken: ", time.strftime("%H:%M:%S",time.gmtime(seconds)), "\n")

