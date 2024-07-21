import os
from stage.stage_builder import StageBuilder


class CutPrimerStage(StageBuilder):
    def __init__(self, config, heading="cut_primer", merge_dir="", save_dir="",
                 rm_p_5="GTCGGTAAAACTCGTGCCAGC",
                 rm_p_3="CAAACTGGGATTAGATACCCCACTATG",
                 error_rate=0.15,
                 min_read_len=204,
                 max_read_len=254,
                 ):
        super().__init__(heading=heading, config=config)
        self.CUTADAPT_PROG = "cutadapt"
        self.merge_dir = merge_dir
        self.save_dir = save_dir
        self.parse_params(rm_p_5, rm_p_3, min_read_len, max_read_len, error_rate)

    def parse_params(self, rm_p_5, rm_p_3, min_read_len, max_read_len, error_rate):
        adapter_length = len(rm_p_5) + len(rm_p_3)
        self.params = (
            f'-g "{rm_p_5};max_error_rate={error_rate}...{rm_p_3};max_error_rate={error_rate}"'
            f" --minimum-length {min_read_len - adapter_length}"
            f" --maximum-length {max_read_len - adapter_length}"
        )

    def setup(self, prefix):
        infile = os.path.join(self.merge_dir, f"{prefix}_merge.fastq")
        cutprimer_outfile = os.path.join(self.save_dir, f"{prefix}_cut.fastq")
        report = os.path.join(self.save_dir, f"{prefix}_report.txt")
        cmd = (
            f"{self.CUTADAPT_PROG} {infile}"
            f" {self.params}"
            f" --discard-untrimmed -j {self.config.n_cpu}"
            f" >{cutprimer_outfile} 2>{report}"
        )
        super().add_stage("cutadapt", cmd)
        # super().add_stage_output_to_file("Redirect cutadapt output", 0, cutprimer_outfile)

    def run(self):
        super().run()
        return all(self.output)


def cutadapt_demo(config, prefix, merge_dir="", save_dir=""):
    stage = CutPrimerStage(config, merge_dir=merge_dir, save_dir=save_dir)
    stage.setup(prefix)
    is_complete = stage.run()
    return is_complete
