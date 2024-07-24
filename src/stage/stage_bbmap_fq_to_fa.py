import os
from edna_processor.utils.base_logger import logger
from stage.stage_builder import StageBuilder
from stage.stage_config import StageConfig

class FqToFaStage(StageBuilder):
    def __init__(self, config, heading="fq_to_fa", cutprimer_dir="", save_dir="",
                 overwrite=True
        ):
        super().__init__(heading=heading, config=config)
        self.BBMAP_PROG = "reformat.sh"
        self.in_suffix = "cut.fastq"
        self.out_suffix = "cut.fasta"
        self.cutprimer_dir = cutprimer_dir
        self.save_dir = save_dir
        self.parse_params(overwrite)

    def parse_params(self, overwrite):
        self.params = (
            f"overwrite={overwrite}"
        )

    def setup(self, prefix):
        infile = os.path.join(self.cutprimer_dir, f"{prefix}_{self.in_suffix}").replace("\\", "/")
        fasta_outfile = os.path.join(self.save_dir, f"{prefix}_{self.out_suffix}").replace("\\", "/")
        report = os.path.join(self.save_dir, f"{prefix}_report.txt").replace("\\", "/")
        self.check_path(infile)
        cmd = (
            f"bash {self.BBMAP_PROG}"
            f" in={infile} out={fasta_outfile}"
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