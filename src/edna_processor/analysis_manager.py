from edna_processor.data_container import SamplesContainer

class AnalysisManager(SamplesContainer):
    def __init__(self, load_path=None):
        super().__init__(load_path)
        self.analysis_type = None
        self.sample_id_used = None
        self.parameters = {}
        self.results_dir = None


    def load_sample_id_list(self, sample_id_list: str = []) -> list[str]:
        if sample_id_list == []:
            sample_id_list = self.sample_id_list
            print(f"> No sample ID list specified. Using all {len(sample_id_list)} samples.")
        else:
            print(f"> Specified {len(sample_id_list)} samples.")
        self.sample_id_used = sample_id_list
        return sample_id_list
    
    def load_units2fasta(self,
            target_name: str,
            target_level: str,
            unit_level: str,
            sample_id_list: list[str]
        ) -> dict[str, str]:
        units2fasta = {}
        for sample_id in sample_id_list:
            for hap, level_dict in self.sample_data[sample_id].hap2level.items():
                if target_name not in level_dict[target_level]:
                    continue
                unit_name = level_dict[unit_level]
                title = f"{unit_name}-{sample_id}_{hap}"
                seq = self.sample_data[sample_id].hap_seq[hap]

                if unit_name not in units2fasta:
                    units2fasta[unit_name] = ""
                units2fasta[unit_name] += f'>{title}\n{seq}\n'
        return units2fasta