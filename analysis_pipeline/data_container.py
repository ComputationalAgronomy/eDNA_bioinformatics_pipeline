import os
import pickle
import time
from datetime import date

from read_blast_csv import BlastReader
from read_denoise_report import DenoiseReportReader
from read_fasta import FastaReader

class OneSampleContainer():
    """
    Container for handling and processing data from various bioinformatics files.

    This class initializes by reading data from given FASTA, denoise report, and BLAST table files,
    and stores the parsed information for further use.

    :attribute amp_seq: A dictionary containing amplicon sequences from the unique FASTA file.
    :attribute hap_seq: A dictionary containing haplotype sequences from the ZOTU FASTA file.
    :attribute amp_size: A dictionary containing amplicon sizes from the denoise report.
    :attribute hap2amp: A dictionary mapping haplotypes to amplicons from the denoise report.
    :attribute hap_size: A dictionary containing haplotype sizes from the denoise report.
    :attribute hap2level: A dictionary mapping haplotypes to taxonomic levels from the BLAST table.

    :param uniq_fasta_path: Path to the unique amplicon FASTA file.
    :param zotu_fasta_path: Path to the ZOTU haplotype FASTA file.
    :param denoise_report_path: Path to the denoise report file.
    :param blast_table_path: Path to the BLAST table file.
    """
    def __init__(self,
            uniq_fasta_path: str,
            zotu_fasta_path: str,
            denoise_report_path: str,
            blast_table_path: str
        ):
        ufr = FastaReader()
        ufr.read_fasta(seq_path=uniq_fasta_path, seq_type="Amplicon")
        self.amp_seq = ufr.seq_dict

        zfr = FastaReader()
        zfr.read_fasta(seq_path=zotu_fasta_path, seq_type="Haplotype")
        self.hap_seq = zfr.seq_dict

        drr = DenoiseReportReader()
        drr.read_denoise_report(denoise_report_path=denoise_report_path)
        self.amp_size = drr.amp_size
        self.hap2amp = drr.hap2amp
        self.hap_size = drr.hap_size

        br = BlastReader()
        br.read_blast_table(blast_table_path=blast_table_path)
        self.hap2level = br.hap2level

class SamplesContainer():

    DATA_FILE_INFO = {
        "uniq_fasta": ("4_derep", "_uniq.fasta"),
        "zotu_fasta": ("5_denoise", "_zotu.fasta"),
        "denoise_report": ("5_denoise", "_zotu_report.csv"),
        "blast_table": ("6_blast", "_blast.csv")
    }

    @staticmethod
    def save_instance(instance, path):
        with open(path, 'wb') as f:
            pickle.dump(instance, f)

    @staticmethod
    def get_sample_id_list(parent_dir):
        child_dir, suffix = tuple(SamplesContainer.DATA_FILE_INFO.values())[0]
        child_dir_path = os.path.join(parent_dir, child_dir)
        print(f"> Searching files with the suffix '{suffix}' from the directory: {child_dir_path}.")
        file_list = os.listdir(child_dir_path)
        sample_id_list = [file.replace(suffix, '') for file in file_list if file.endswith(suffix)]
        print(f"> Found {len(sample_id_list)} samples.")
        return sample_id_list

    @staticmethod
    def check_dir(parent_dir):
        for child_dir, _ in SamplesContainer.DATA_FILE_INFO.values():
            child_dir_path = os.path.join(parent_dir, child_dir)
            if not os.path.isdir(child_dir_path):
                print(f"> Error: Directory {child_dir_path} does not exist")
                return 
        print("> All directories exist.")

    @staticmethod
    def check_file(parent_dir, sample_id_list):
        print("\n> Start checking whether files exist...")
        for child_dir, suffix in SamplesContainer.DATA_FILE_INFO.values():
            for sample_id in sample_id_list:
                file_path = os.path.join(parent_dir, child_dir, f"{sample_id}{suffix[0]}")
                if not os.path.isfile(file_path):
                    print(f"> Error: {file_path} does not exist.")
                    return
        print("> All files exist.")

    @staticmethod
    def get_file_paths(sample_id, parent_dir):
        file_paths = {}
        for file_key, (child_dir, suffix) in SamplesContainer.DATA_FILE_INFO.items():
            file_path = os.path.join(parent_dir, child_dir, f"{sample_id}{suffix}")
            file_paths[file_key] = file_path
        return file_paths

    def __init__(self, load_path=None):
        self.sample_data = {}
        self.sample_id_list = []

        if load_path:
            self.load_sample_data(load_path)

    def import_samples(self, parent_dir, sample_id_list=None): # read_dir: path to the directory containing the "4_derep", "5_denoise", "6_blast" folders.

        start_time = time.time()

        print(f"> Reading samples from: {parent_dir}.")
        SamplesContainer.check_dir(parent_dir)

        if sample_id_list is None:
            print(f"> No sample id list provided.")
            SamplesContainer.get_sample_id_list(parent_dir)
            SamplesContainer.check_file(parent_dir, sample_id_list)
        else:
            print("> Specified sample id list.")
            SamplesContainer.check_file(parent_dir, sample_id_list)

        self.sample_id_list.extend(sample_id_list)

        for sample_id in sample_id_list:
            file_paths = SamplesContainer.get_file_paths(sample_id, parent_dir)

            uniq_fasta_path = file_paths["uniq_fasta"]
            zotu_fasta_path = file_paths["zotu_fasta"]
            zotu_report_path = file_paths["denoise_report"]
            blast_table_path = file_paths["blast_table"]

            self.sample_data[sample_id] = OneSampleContainer(uniq_fasta_path=uniq_fasta_path, zotu_fasta_path=zotu_fasta_path, denoise_report_path=zotu_report_path, blast_table_path=blast_table_path)

        seconds = time.time() - start_time
        print(f"> {len(sample_id_list)} Samples read.")
        print("Time Taken: ", time.strftime("%H:%M:%S",time.gmtime(seconds)), "\n")

    def save_sample_data(self, save_dir, save_name=f"eDNA_samples_{date.today()}"):
        os.makedirs(save_dir, exist_ok=True)

        save_path = os.path.join(save_dir, f"{save_name}.pkl")
        print(f"> Saving sample data to: {save_path}.")

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

        SamplesContainer.save_instance(self, save_path)
        print(f"> Sample data saved to: {save_path}\n")

    def load_sample_data(self, load_path):
        with open(load_path,'rb') as file:
            self.__dict__ = pickle.load(file).__dict__
        print(f"> Sample data loaded from: {load_path}\n")