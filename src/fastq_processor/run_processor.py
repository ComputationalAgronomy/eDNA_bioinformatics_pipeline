import os
import sys

from analysis_toolkit.utils.base_logger import logger
from fastq_processor.step_exec.stage_bbmap_fq_to_fa import FqToFaStage
from fastq_processor.step_exec.stage_blastn_assign_taxa import AssignTaxaStage
from fastq_processor.step_build.stage_config import StageConfig
from fastq_processor.step_exec.stage_cutadapt_cut_primer import CutPrimerStage
from fastq_processor.step_exec.stage_usearch_dereplicate import DereplicateStage
from fastq_processor.step_exec.stage_usearch_denoise import DenoiseStage
from fastq_processor.step_exec.stage_usearch_merge import MergeStage


class FastqProcessor:
    
    @staticmethod
    def get_prefix_with_suffix(in_dir, suffix):
        files = os.listdir(in_dir)
        prefix = [file.replace(suffix, "") for file in files if file.endswith(suffix)]
        return prefix

    @staticmethod
    def setup_stages(config, stages_parent_dir, fastq_dir_name, merge_dir_name, cutprimer_dir_name, fqtofa_dir_name, derep_dir_name, denoise_dir_name, blast_dir_name, db_path, lineage_path):
        fastq_dir = os.path.join(stages_parent_dir, fastq_dir_name)
        merge_dir = os.path.join(stages_parent_dir, merge_dir_name)
        cutprimer_dir = os.path.join(stages_parent_dir, cutprimer_dir_name)
        fqtofa_dir = os.path.join(stages_parent_dir, fqtofa_dir_name)
        derep_dir = os.path.join(stages_parent_dir, derep_dir_name)
        denoise_dir = os.path.join(stages_parent_dir, denoise_dir_name)
        blast_dir = os.path.join(stages_parent_dir, blast_dir_name)

        stages = dict()
        stages["merge"] = MergeStage(config, fastq_dir=fastq_dir, save_dir=merge_dir)
        stages["cutprimer"] = CutPrimerStage(config, merge_dir=merge_dir, save_dir=cutprimer_dir)
        stages["fqtofa"] = FqToFaStage(config, cutprimer_dir=cutprimer_dir, save_dir=fqtofa_dir)
        stages["dereplicate"] = DereplicateStage(config, fasta_dir=fqtofa_dir, save_dir=derep_dir)
        stages["denoise"] = DenoiseStage(config, derep_dir=derep_dir, save_dir=denoise_dir)
        stages["assigntaxa"] = AssignTaxaStage(config, denoise_dir=denoise_dir, save_dir=blast_dir,
                                               db_path=db_path, lineage_path=lineage_path)

        return stages

    @staticmethod
    def run_each_data(prefix, stages):
        for k, s in stages.items():
            print(prefix)
            s.setup(prefix)
            # print(s.runners[0].command)
            is_complete = s.run()
            if not is_complete:
                print(f"Error: process errors at stage: {k}\n")
                break

    def __init__(self,
        stages_parent_dir: str,
        fastq_dir_name: str,
        db_path: str,
        lineage_path: str,
        merge_dir_name: str = "merge",
        cutprimer_dir_name: str = "cut_primer",
        fqtofa_dir_name: str = "fq_to_fa",
        derep_dir_name: str = "dereplicate",
        denoise_dir_name: str = "denoise",
        blast_dir_name: str = "blast",
        verbose: bool=True,
        n_cpu: int = 1,
        memory: int = 8,
        ):
        self.config = StageConfig(verbose=verbose, logger=logger, n_cpu=n_cpu, memory=memory)
        self.parent_dir = stages_parent_dir
        self.input_dir = os.path.join(stages_parent_dir, fastq_dir_name)
        self.stages = FastqProcessor.setup_stages(
            self.config,
            stages_parent_dir,
            fastq_dir_name,
            merge_dir_name,
            cutprimer_dir_name,
            fqtofa_dir_name,
            derep_dir_name,
            denoise_dir_name,
            blast_dir_name,
            db_path,
            lineage_path,
        )
        self.data_prefix = FastqProcessor.get_prefix_with_suffix(self.input_dir, "_R1.fastq")
        for prefix in self.data_prefix:
            FastqProcessor.run_each_data(prefix, self.stages)

def main():
    FastqProcessor(
        stages_parent_dir="stage_test",
        fastq_dir_name="fastq",
        db_path="data\\database\\MiFish",
        lineage_path="data\\database\\lineage.csv",
        n_cpu=20,
    )

if __name__ == "__main__":
    main()