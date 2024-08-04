import os

from fastq_processor.step_build.stage_builder import StageBuilder


class MergeStage(StageBuilder):
    def __init__(self, config, heading="stage_usearch_merge.py", fastq_dir="", save_dir="",
                 maxdiff=5,
                 pctid=90
        ):
        super().__init__(heading=heading, config=config)
        self.USEARCH_PROG = "usearch" # TODO(SW): Don't use `.exe`, doesn't make sense in docker/ubuntu
        self.in_suffix = "R1.fastq"
        self.out_suffix = "merge.fastq"
        self.report_suffix = "report.txt"
        self.fastq_dir = fastq_dir
        self.save_dir = save_dir
        self.parse_params(maxdiff, pctid)

    def parse_params(self, maxdiff, pctid):
        self.params = (
            f"-fastq_maxdiffs {maxdiff} -fastq_pctid {pctid}"
            f" -threads {self.config.n_cpu}"
        )

    def setup(self, prefix):
        infile = os.path.join(self.fastq_dir, f"{prefix}_{self.in_suffix}")
        merge_outfile = os.path.join(self.save_dir, f"{prefix}_{self.out_suffix}")
        report = os.path.join(self.save_dir, f"{prefix}_{self.report_suffix}")
        self.check_path(infile)
        cmd = (
            f"{self.USEARCH_PROG}"
            f" -fastq_mergepairs {infile} -fastqout {merge_outfile}"
            f" {self.params}"
            f" -report {report}"
        )
        super().add_stage("Merge paired-end sequences", cmd)

    def run(self):
        super().run()
        return all(self.output)


def usearch_merge_demo(config, prefix, fastq_dir, save_dir):
    stage = MergeStage(config, fastq_dir=fastq_dir, save_dir=save_dir)
    stage.setup(prefix)
    is_complete = stage.run()
    return is_complete