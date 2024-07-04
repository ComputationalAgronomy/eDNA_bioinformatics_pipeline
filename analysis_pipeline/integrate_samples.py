import os
import pickle
import time
from datetime import date

from read_blast_csv import read_hap_level
from read_denoise_report import read_denoise_report
from read_fasta import read_fasta_file


def save_instance(instance, path):
    with open(path, 'wb') as f:
        pickle.dump(instance, f)

def get_sample_id_list(read_dir):
    file_list = os.listdir(f'{read_dir}/4_derep')
    sample_id_list = [file.replace('_uniq.fasta', '') for file in file_list if file.endswith('_uniq.fasta')]
    return sample_id_list

def check_dir(read_dir):
    # TODO(SW): simplify this with ["4_derep", "5_denoise", "6_blast"]
    if os.path.exists(os.path.join(read_dir, "4_derep")) is False:
        print("> Error: 4_derep folder does not exist in the specified directory.")
        return
    if os.path.exists(os.path.join(read_dir, "5_denoise")) is False:
        print("> Error: 5_denoise folder does not exist in the specified directory.")
        return
    if os.path.exists(os.path.join(read_dir, "6_blast")) is False:
        print("> Error: 6_blast folder does not exist in the specified directory.")
        return

def check_file(read_dir, sample_id_list):
    for sample_id in sample_id_list:
        if os.path.exists(os.path.join(read_dir, "4_derep", f"{sample_id}_uniq.fasta")) is False:
            print(f"> Error: {sample_id}_uniq.fasta does not exist in the specified directory.")
            return
        if os.path.exists(os.path.join(read_dir, "5_denoise", f"{sample_id}_zotu.fasta")) is False:
            print(f"> Error: {sample_id}_zotu.fasta does not exist in the specified directory.")
            return
        if os.path.exists(os.path.join(read_dir, "6_blast", f"{sample_id}_blast.csv")) is False:
            print(f"> Error: {sample_id}_blast.csv does not exist in the specified directory.")
            return
    print("> All files exist.")

class OneSample():

    def __init__(self, uniq_fasta_path, zotu_fasta_path, zotu_report_path, blast_csv_path):
        self.amp_seq = read_fasta_file(uniq_fasta_path, seq_type="Amplicon")
        self.amp_size, self.hap2amp, self.hap_size = read_denoise_report(zotu_report_path)
        self.hap_seq = read_fasta_file(zotu_fasta_path, seq_type="Haplotype")
        self.hap2level = read_hap_level(blast_csv_path)

class IntegrateSamples():

    def __init__(self, load_path=None):
        self.sample_data = {}
        self.sample_id_list = []

        if load_path:
            self.load_sample_data(load_path)

    def import_samples(self, read_dir, sample_id_list=None): # read_dir: path to the directory containing the "4_derep", "5_denoise", "6_blast" folders.

        start_time = time.time()

        print(f"> Reading samples from: {read_dir}.")
        check_dir(read_dir)

        if sample_id_list is None:
            print("> No sample id list provided, reading sample id list from the directory.")
            sample_id_list = get_sample_id_list(read_dir)
        else:
            print("> Specified sample id list. Start checking whether files exist.")
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

    def save_sample_data(self, save_dir, save_name=f"eDNA_samples_{date.today()}"):
        print(f"> Saving sample data to: {save_dir}.")

        save_path = os.path.join(save_dir, f"{save_name}.pkl")
        # overwrite_y_n
        if os.path.exists(save_path):
            while True:
                user_input = input("> File already exists, Do you want to overwrite it? (y/n)")
                if user_input in ['y', 'Y', 'Yes', 'yes']:
                    os.remove(save_path)
                    print("> File overwritten.")
                    break
                elif user_input in ['n', 'N', 'No', 'no']:
                    print("> File not saved.\n")
                    return
                else:
                    print("> Invalid input.")

        save_instance(self, save_path)
        print(f"> Sample data saved to: {save_path}\n")

    def load_sample_data(self, load_path):
        with open(load_path,'rb') as file:
            self.__dict__ = pickle.load(file).__dict__
        print(f"> Sample data loaded from: {load_path}\n")