import os

from fastq_processor.step_build.stage_builder import StageBuilder


class FqToFaStage(StageBuilder):
    def __init__(self, config, heading="stage_bbmap_fq_to_fa.py", cutprimer_dir="", save_dir="",
                 overwrite=True
        ):
        super().__init__(heading=heading, config=config)
        self.BBMAP_PROG = "reformat.sh"
        self.in_suffix = "cut.fastq"
        self.out_suffix = "cut.fasta"
        self.report_suffix = "report.txt"
        self.cutprimer_dir = cutprimer_dir
        self.save_dir = save_dir
        self.parse_params(overwrite)

    def parse_params(self, overwrite):
        self.params = (
            f"overwrite={overwrite}"
        )

    def setup(self, prefix):
        self.infile = os.path.join(self.cutprimer_dir, f"{prefix}_{self.in_suffix}").replace("\\", "/")
        fasta_outfile = os.path.join(self.save_dir, f"{prefix}_{self.out_suffix}").replace("\\", "/")
        report = os.path.join(self.save_dir, f"{prefix}_{self.report_suffix}")
        self.check_infile()
        self.check_savedir()
        cmd = (
            f"bash {self.BBMAP_PROG}"
            f" in={self.infile} out={fasta_outfile}"
            f" {self.params}"
        )
        super().add_stage("Reformat FASTQ to FASTA", cmd)
        super().add_stage_output_to_file("Write reformat.sh report", 0, report, report)

    def run(self):
        super().run()
        return all(self.output)


def fq_to_fa_demo(config, prefix, cutprimer_dir="", save_dir=""):
    stage = FqToFaStage(config, cutprimer_dir=cutprimer_dir, save_dir=save_dir)
    stage.setup(prefix)
    is_complete = stage.run()
    return is_complete