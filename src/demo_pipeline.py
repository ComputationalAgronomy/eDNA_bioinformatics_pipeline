import os
import sys
from stage.cutadapt_stage import CutadaptStage
from stage.usearch_merge_stage import UsearchMergeStage
from stage.stage_config import StageConfig


def get_prefix_with_suffix(in_dir, suffix):
    files = os.listdir(in_dir)
    prefix = []
    for filename in files:
        if filename.endswith(suffix):
            prefix.append(filename.replace(suffix, ""))
    return prefix


def setup_stages(config, in_dir):
    fastq_dir = in_dir
    merge_dir = "merge_dir"
    save_dir = "output"
    stages = dict()
    stages["usearch"] = UsearchMergeStage(config, fastq_dir=fastq_dir, save_dir=save_dir)
    stages["cutadapt"] = CutadaptStage(config, merge_dir=merge_dir, save_dir=save_dir)
    # stages["fq_to_fa"] = FqToFaStage(config, save_dir=save_dir)
    # stages["dereplicate"] = UsearchDereplicateStage(config, save_dir=save_dir)
    return stages


def run_each_data(prefix, stages):
    for k, s in stages.items():
        print(prefix)
        s.setup(prefix)
        print(s.runners[0].command)
        is_complete = s.run()
        if not is_complete:
            print(f"Error: process errors at stage: {k}\n")
            break


def main():
    in_dir = "data"
    config = StageConfig(verbose=True, dry=True, logger=sys.stdout, n_cpu=1, memory=8)
    stages = setup_stages(config, in_dir)
    data_prefix = get_prefix_with_suffix(in_dir, "_R1.fastq")
    for prefix in data_prefix:
        run_each_data(prefix, stages)

if __name__ == "__main__":
    main()