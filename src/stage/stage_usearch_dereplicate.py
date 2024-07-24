import os
from stage.stage_builder import StageBuilder
from edna_processor.utils.base_logger import logger
from stage.stage_config import StageConfig


class DereplicateStage(StageBuilder):
    def __init__(self, config, heading="dereplicate", fasta_dir="", save_dir="",
                 annot_size: bool = True,
                 seq_label: str = "Uniq"
        ):
        super().__init__(heading=heading, config=config)
        self.USEARCH_PROG = "usearch.exe"
        self.in_suffix = "cut.fasta"
        self.out_suffix = "derep.fasta"
        self.fasta_dir = fasta_dir
        self.save_dir = save_dir
        self.parse_params(annot_size, seq_label)

    def parse_params(self, annot_size, seq_label):
        sizeout = "-sizeout" if annot_size else ""
        self.params = (
            f"{sizeout} -relabel {seq_label}"
        )

    def setup(self, prefix):
        infile = os.path.join(self.fasta_dir, f"{prefix}_{self.in_suffix}")
        dereplicate_outfile = os.path.join(self.save_dir, f"{prefix}_{self.out_suffix}")
        report = os.path.join(self.save_dir, f"{prefix}_report.txt")
        self.check_path(infile)
        cmd = (
            f"{self.USEARCH_PROG} -fastx_uniques {infile}"
            f" {self.params}"
            f" -threads {self.config.n_cpu}"
            f" -fastaout {dereplicate_outfile}"
        )
        super().add_stage("Dereplicate trimmed sequences", cmd)
        super().add_stage_output_to_file("Write usearch report", 0, report, report)

    def run(self):
        super().run()
        return all(self.output)


def dereplicate_demo(config, prefix, fasta_dir="", save_dir=""):
    stage = DereplicateStage(config, fasta_dir=fasta_dir, save_dir=save_dir)
    stage.setup(prefix)
    is_complete = stage.run()
    return is_complete