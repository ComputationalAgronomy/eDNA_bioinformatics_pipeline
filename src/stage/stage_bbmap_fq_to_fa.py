import os
from stage.stage_builder import StageBuilder


class FqToFaStage(StageBuilder):
    def __init__(self, config, heading="fq_to_fa", cutprimer_dir="", save_dir="",
                 overwrite=True
        ):
        super().__init__(heading=heading, config=config)
        self.BBMAP_PROG = "reformat.sh"
        self.cutprimer_dir = cutprimer_dir
        self.save_dir = save_dir
        self.parse_params(overwrite)

    def parse_params(self, overwrite):
        self.params = (
            f"overwrite={overwrite}"
        )

    def setup(self, prefix):
        infile = os.path.join(self.cutprimer_dir, f"{prefix}_cut.fastq")
        fasta_outfile = os.path.join(self.save_dir, f"{prefix}_cut.fasta")
        report = os.path.join(self.save_dir, f"{prefix}_report.txt")
        cmd = (
            f"bash {self.BBMAP_PROG}"
            f" in={infile} out={fasta_outfile}"
            f" {self.params}"
            f" 2>{report}"
        )
        super().add_stage("bbmap_reformat.sh", cmd)

    def run(self):
        super().run()
        return all(self.output)


def fq_to_fa_demo(config, prefix, cutprimer_dir="", save_dir=""):
    stage = FqToFaStage(config, cutprimer_dir=cutprimer_dir, save_dir=save_dir)
    stage.setup(prefix)
    is_complete = stage.run()
    return is_complete