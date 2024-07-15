import os
import pickle
from datetime import date

from edna_processor.read.read_blast_csv import BlastReader
from edna_processor.read.read_denoise_report import DenoiseReportReader
from edna_processor.read.read_fasta import FastaReader
from edna_processor.utils.base_logger import logger

class OneSampleContainer():
    """
    Container for handling and processing data from various bioinformatics files.
    This class initializes by reading data from given FASTA, denoise report, and BLAST table files,
    and stores the parsed information for further use.

    :param uniq_fasta_path: Path to the unique amplicon FASTA file.
    :param zotu_fasta_path: Path to the ZOTU haplotype FASTA file.
    :param denoise_report_path: Path to the denoise report file.
    :param blast_table_path: Path to the BLAST table file.

    :attribute amp_seq: A dictionary containing amplicon sequences from the unique FASTA file.
    :attribute hap_seq: A dictionary containing haplotype sequences from the ZOTU FASTA file.
    :attribute amp_size: A dictionary containing amplicon sizes from the denoise report.
    :attribute hap2amp: A dictionary mapping haplotypes to amplicons from the denoise report.
    :attribute hap_size: A dictionary containing haplotype sizes from the denoise report.
    :attribute hap2level: A dictionary mapping haplotypes to taxonomic levels from the BLAST table.

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
    """
    A class to manage sample data files.

    :attribute DATA_FILE_INFO: A dictionary mapping input data types to tuples containing the corresponding child directory and file suffix.
    :attribute sample_data: A dictionary to store sample data.
    :attribute sample_id_list: A list to store sample IDs.
    :attribute data_parent_dir: The parent directory path where the data files are located.
    :attribute instance_path: The file path where the instance of SamplesContainer is saved or loaded.
    """
    # 
    # keys are data type, values are tuples of (child directory, file suffix).
    DATA_FILE_INFO = {
        "uniq_fasta": ("4_derep", "_uniq.fasta"),
        "zotu_fasta": ("5_denoise", "_zotu.fasta"),
        "denoise_report": ("5_denoise", "_zotu_report.csv"),
        "blast_table": ("6_blast", "_blast.csv")
    }

    @staticmethod
    def save_instance(instance: 'SamplesContainer', path: str) -> None:
        """
        Save an instance of SamplesContainer to a specified path.

        :param instance: The SamplesContainer instance to save.
        :param path: The file path where the instance will be saved.
        """
        with open(path, 'wb') as f:
            pickle.dump(instance, f)

    @staticmethod
    def get_sample_id_list(parent_dir: str) -> list[str]:
        """
        Retrieve a list of sample IDs from the specified parent directory.

        :param parent_dir: Path to the parent directory containing sample data.
        :return: A list of sample IDs.
        """
        child_dir, suffix = tuple(SamplesContainer.DATA_FILE_INFO.values())[0]

        child_dir_path = os.path.join(parent_dir, child_dir)

        logger.info(f"Searching files with the suffix '{suffix}' from the directory: {child_dir_path}.")
    
        file_list = os.listdir(child_dir_path)
        sample_id_list = [file.replace(suffix, '') for file in file_list if file.endswith(suffix)]

        logger.info(f"Found {len(sample_id_list)} samples.")
    
        return sample_id_list

    @staticmethod
    def check_dir(parent_dir: str) -> None:
        """
        Check if all necessary child directories exist within the specified parent directory.

        :param parent_dir: Path to the parent directory to check.
        """
        for child_dir, _ in SamplesContainer.DATA_FILE_INFO.values():
            child_dir_path = os.path.join(parent_dir, child_dir)

            if not os.path.isdir(child_dir_path):
                raise FileNotFoundError(f"Directory does not exist: {child_dir_path}.")

        logger.info("All directories exist.")

    @staticmethod
    def check_file(parent_dir: str, sample_id_list: list[str]) -> None:
        """
        Check if all necessary files for the given sample IDs exist within the specified parent directory.

        :param parent_dir: Path to the parent directory to check.
        :param sample_id_list: List of sample IDs to check.
        """
        logger.info("Start checking whether files exist...")

        for child_dir, suffix in SamplesContainer.DATA_FILE_INFO.values():
            for sample_id in sample_id_list:
                file_path = os.path.join(parent_dir, child_dir, f"{sample_id}{suffix[0]}")

                if not os.path.isfile(file_path):
                    raise FileNotFoundError(f"File does not exist: {file_path}.")

        logger.info("All files exist.")

    @staticmethod
    def get_file_paths(sample_id: str, parent_dir: str) -> dict[str, str]:
        """
        Retrieve the file paths for all data types for a given sample ID.

        :param sample_id: The sample ID to retrieve file paths for.
        :param parent_dir: Path to the parent directory containing sample data.
        :return: A dictionary mapping data types to their respective file paths.
        """
        file_paths = {}
        for file_key, (child_dir, suffix) in SamplesContainer.DATA_FILE_INFO.items():
            file_path = os.path.join(parent_dir, child_dir, f"{sample_id}{suffix}")
            file_paths[file_key] = file_path

        return file_paths

    def __init__(self, load_path: str = None):
        """
        Initialize a SamplesContainer instance.

        :param load_path: If provided, load a pre-existing SamplesContainer instance from the path. Defaults is None.
        """
        self.sample_data = {}
        self.sample_id_list = []
        self.data_imported_dir = None
        self.instance_path = None

        if load_path:
            self.load_sample_data(load_path)

    def import_samples(self, parent_dir: str, sample_id_list: list[str] = None) -> None:
        """
        Import sample data from the specified parent directory.
        The parent directory is expected to contain four data types recorded in 'DATA_FILE_INFO'.
        Data should be organized in child directories with specific suffixes as defined in 'DATA_FILE_INFO'.

        :param parent_dir: Path to the parent directory containing the sample data.
        :param sample_id_list: List of sample IDs to import. If not provided, all available sample IDs will be imported. The sample IDs are extracted from the file names using the provided suffix.
        """
        logger.info(f"Reading samples from: {parent_dir}.")
        SamplesContainer.check_dir(parent_dir)

        if sample_id_list is None:
            logger.info(f"No sample id list provided.")
            SamplesContainer.get_sample_id_list(parent_dir)
            SamplesContainer.check_file(parent_dir, sample_id_list)
        else:
            logger.info("Specified sample id list.")
            SamplesContainer.check_file(parent_dir, sample_id_list)

        self.data_imported_dir = parent_dir
        self.sample_id_list.extend(sample_id_list)

        for sample_id in sample_id_list:
            file_paths = SamplesContainer.get_file_paths(sample_id, parent_dir)

            uniq_fasta_path = file_paths["uniq_fasta"]
            zotu_fasta_path = file_paths["zotu_fasta"]
            zotu_report_path = file_paths["denoise_report"]
            blast_table_path = file_paths["blast_table"]

            self.sample_data[sample_id] = OneSampleContainer(uniq_fasta_path=uniq_fasta_path, zotu_fasta_path=zotu_fasta_path, denoise_report_path=zotu_report_path, blast_table_path=blast_table_path)

        logger.info(f"{len(sample_id_list)} Samples read.")

    def save_sample_data(self, save_dir: str, save_name: str = f"eDNA_samples_{date.today()}") -> None:
        """
        Save the current sample data to a specified directory.

        :param save_dir: The directory where the sample data will be saved.
        :param save_name: The name of the save file. Defaults to 'eDNA_samples_<current_date>'.
        """
        os.makedirs(save_dir, exist_ok=True)

        save_path = os.path.join(save_dir, f"{save_name}.pkl")
        logger.info(f"Saving sample data to: {save_path}.")

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
        self.instance_path = save_path
        logger.info(f"Sample data saved to: {save_path}\n")

    def load_sample_data(self, load_path: str) -> None:
        """
        Load sample data from a specified path.

        :param load_path: The file path to load the sample data from.
        """
        with open(load_path,'rb') as file:
            self.__dict__ = pickle.load(file).__dict__
        self.instance_path = load_path
        logger.info(f"Sample data loaded from: {load_path}\n")