from abc import ABC, abstractmethod

from analysis_toolkit import data_container
from analysis_toolkit.utils import base_logger


class Runner(ABC):

    def __init__(self, samplesdata: data_container.SamplesData):
        self.analysis_type = None
        self.sample_id_used = None
        self.parameters = {}
        self.results_dir = None
        self.import_data(samplesdata)

    def import_data(self, samplesdata):
        self.sample_data = samplesdata.sample_data
        self.sample_id_list = samplesdata.sample_id_list

    def load_sample_id_list(self, sample_id_list: str = []):
        if sample_id_list == []:
            base_logger.logger.info(f"No sample ID list specified. Using all {len(sample_id_list)} samples.")
            self.sample_id_used = self.sample_id_list
        else:
            base_logger.logger.info(f"Specified {len(sample_id_list)} samples.")
            for sample_id in sample_id_list:
                if sample_id not in self.sample_id_list:
                    raise ValueError(f"Specified invalid sample ID: {sample_id}.")
            self.sample_id_used = sample_id_list

    # TODO(SW): Refactor **_target() and **write_output() with the following
    @abstractmethod
    def target(self, *args):
        pass

    @abstractmethod
    def write_output(self, *args):
        pass

class SequenceRunner(ABC, Runner):
    def __init__(self, samplesdata: data_container.SamplesData):
        super().__init__(samplesdata)
        self.units2fasta = {}

    def load_units2fasta_dict(self,
            target_name: str,
            target_level: str,
            unit_level: str,
            n_unit_threshold: int = -1 # TODO(SW): OR *args
        ):
        for sample_id in self.sample_id_used:
            for hap, level_dict in self.sample_data[sample_id].hap2level.items():
                if target_name not in level_dict[target_level]:
                    continue

                unit_name = level_dict[unit_level]
                title = f"{unit_name}-{sample_id}_{hap}"
                seq = self.sample_data[sample_id].hap_seq[hap]

                if unit_name not in self.units2fasta:
                    self.units2fasta[unit_name] = ""
                self.units2fasta[unit_name] += f'>{title}\n{seq}\n'

class AbundanceRunner(ABC, Runner):
    def __init__(self, samplesdata: data_container.SamplesData):
        super().__init__(samplesdata)
        self.units2abundance = {}

    def load_hap_size(self, sample_id: str, hap: str) -> int:
        """
        Get the size of a haplotype for a given sample.

        :param sample_id: The ID of the sample.
        :param hap: The ID of the haplotype.
        :return: The size of the haplotype.
        """
        return int(self.sample_data[sample_id].hap_size[hap])

    def load_units2abundance_dict(self, sample_id: str, unit_level: str):
        """
        Get the abundance of a level for a given sample.

        :param sample_id: The ID of the sample.
        :param level: The name of the level.
        :return: A dictionary mapping level names to abundances.
        """
        for hap, level_dict in self.sample_data[sample_id].hap2level.items():
            level_name = level_dict[unit_level]
            if level_name not in self.units2abundance:
                self.units2abundance[level_name] = 0
            size = self.load_hap_size(sample_id, hap)
            self.units2abundance[level_name] += size # e.g. {'SpA': 3, 'SpB': 4, 'SpC': 5}
    
    def normalize_abundance(self):
        """
        Normalize the abundance values in the dictionary to percentages.

        :param abundance_dict: A dictionary with unit names and their abundance.
        :returns: A dictionary with unit names and their normalized abundance in percentages.
        """
        total_size = sum(self.units2abundance.values())
        self.units2abundance = {key: value/total_size * 100 for key, value in self.units2abundance.items()}
