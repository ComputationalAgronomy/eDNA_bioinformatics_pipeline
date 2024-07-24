import os
import sys

from edna_processor.utils.base_logger import logger
from stage.stage_bbmap_fq_to_fa import FqToFaStage
from stage.stage_blastn_assign_taxa import AssignTaxaStage
from stage.stage_config import StageConfig
from stage.stage_cutadapt_cut_primer import CutPrimerStage
from stage.stage_usearch_dereplicate import DereplicateStage
from stage.stage_usearch_denoise import DenoiseStage
from stage.stage_usearch_merge import MergeStage


def get_prefix_with_suffix(in_dir, suffix):
    files = os.listdir(in_dir)
    prefix = [file.replace(suffix, "") for file in files if file.endswith(suffix)]
    return prefix


def setup_stages(config, parent_dir):
    fastq_dir = os.path.join(parent_dir, "fastq")
    merge_dir = os.path.join(parent_dir, "merge")
    cutprimer_dir = os.path.join(parent_dir, "cutprimer")
    fqtofa_dir = os.path.join(parent_dir, "fqtofa")
    derep_dir = os.path.join(parent_dir, "dereplicate")
    denoise_dir = os.path.join(parent_dir, "denoise")
    blast_dir = os.path.join(parent_dir, "blast")

    db_path = os.path.join("data", "database", "MiFish")
    lineage_path = os.path.join("data", "database", "lineage.csv")

    stages = dict()
    stages["merge"] = MergeStage(config, fastq_dir=fastq_dir, save_dir=merge_dir)
    stages["cutprimer"] = CutPrimerStage(config, merge_dir=merge_dir, save_dir=cutprimer_dir)
    stages["fqtofa"] = FqToFaStage(config, cutprimer_dir=cutprimer_dir, save_dir=fqtofa_dir)
    stages["dereplicate"] = DereplicateStage(config, fasta_dir=fqtofa_dir, save_dir=derep_dir)
    stages["denoise"] = DenoiseStage(config, derep_dir=derep_dir, save_dir=denoise_dir)
    stages["assigntaxa"] = AssignTaxaStage(config, denoise_dir=denoise_dir, save_dir=blast_dir,
                                           db_path=db_path, lineage_path=lineage_path)

    return stages


def run_each_data(prefix, stages):
    for k, s in stages.items():
        print(prefix)
        s.setup(prefix)
        # print(s.runners[0].command)
        is_complete = s.run()
        if not is_complete:
            print(f"Error: process errors at stage: {k}\n")
            break


def main():
    parent_dir = "stage_test"
    config = StageConfig(verbose=True, dry=False, logger=logger, n_cpu=20, memory=8)
    stages = setup_stages(config, parent_dir)
    data_prefix = get_prefix_with_suffix("stage_test\\fastq","_R1.fastq")
    for prefix in data_prefix:
        run_each_data(prefix, stages)
    # run_each_data("test", stages)

if __name__ == "__main__":
    main()