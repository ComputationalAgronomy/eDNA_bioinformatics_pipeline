import os
from stage.stage_builder import StageBuilder


class DenoiseStage(StageBuilder):
    def __init__(self, config, heading="denoise", derep_dir="", save_dir="",
                 minsize=8,
                 alpha=2
        ):
        super().__init__(heading=heading, config=config)
        self.USEARCH_PROG = "usearch.exe"
        self.derep_dir = derep_dir
        self.save_dir = save_dir
        self.parse_params(minsize, alpha)

    def parse_params(self, minsize, alpha):
        self.params = (
            f"-minsize {minsize} -unoise_alpha {alpha}"
        )

    def setup(self, prefix):
        infile = os.path.join(self.derep_dir, f"{prefix}_derep.fasta")
        denoise_outfile = os.path.join(self.save_dir, f"{prefix}_denoise.fasta")
        denoise_report = os.path.join(self.save_dir, f"{prefix}_denoise_report.txt")
        report = os.path.join(self.save_dir, f"{prefix}_report.txt")
        cmd = (
            f"{self.USEARCH_PROG} -unoise3 {infile}"
            f" {self.params}"
            f" -threads {self.config.n_cpu}"
            f" -zotus {denoise_outfile} -tabbedout {denoise_report}"
            f" >{report} 2>&1"
        )
        super().add_stage("usearch.exe", cmd)

    def run(self):
        super().run()
        return all(self.output)


def fq_to_fa_demo(config, prefix, derep_dir="", save_dir=""):
    stage = DenoiseStage(config, derep_dir=derep_dir, save_dir=save_dir)
    stage.setup(prefix)
    is_complete = stage.run()
    return is_complete