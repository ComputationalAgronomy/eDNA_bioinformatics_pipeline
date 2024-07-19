import os

from stage.stage_builder import StageBuilder


class UsearchMergeStage(StageBuilder):
    def __init__(self, config, heading="usearch",
                 fastq_dir="", save_dir="",
                 maxdiff=5, pctid=90):
        super().__init__(heading=heading, config=config)
        self.USEARCH_PROG = "usearch"
        self.suffix = "R1.fastq"
        self.fastq_dir = fastq_dir
        self.save_dir = save_dir
        self.parse_params(maxdiff, pctid)

    def parse_params(self, maxdiff, pctid):
        self.params = (
            f" -fastq_maxdiffs {maxdiff} -fastq_pctid {pctid}"
            f" -threads {self.config.n_cpu}"
        )

    def setup(self, prefix):

        mergepair = os.path.join(self.fastq_dir, f"{prefix}_{self.suffix}")
        fastqout = os.path.join(self.fastq_dir, f"{prefix}_merged.fastq")
        report = os.path.join(self.save_dir, f"{prefix}_report.txt")
        cmd = (
            f"{self.USEARCH_PROG} "
            f" -fastq_mergepairs {mergepair} -fastqout {fastqout}"
            f"{self.params}"
            f" -report {report}"
        )
        super().add_stage("usearch_merge", cmd, shell=True)

    def run(self):
        super().run()
        return all(self.output)


def usearch_merge_demo(config, prefix, fastq_dir, save_dir):
    stage = UsearchMergeStage(config, fastq_dir=fastq_dir, save_dir=save_dir)
    stage.setup(prefix)
    is_complete = stage.run()
    return is_complete

