import os

from fastq_processor.step_build.stage_builder import StageBuilder


class DecompressStage(StageBuilder):
    def __init__(self, config, heading="stage_gzip_decompress_fastq_gz.py", fastq_dir="", save_dir="",
        ):
        super().__init__(heading=heading, config=config)
        self.BBMAP_PROG = "gzip"
        self.in1_suffix = "R1.fastq.gz"
        self.in2_suffix = "R2.fastq.gz"
        self.out1_suffix = "R1.fastq"
        self.out2_suffix = "R2.fastq"
        self.fastq_dir = fastq_dir
        self.save_dir = save_dir
        self.parse_params()

    def parse_params(self):
        self.params = (
            f"-cd"
        )

    def setup(self, prefix):
        self.infile = os.path.join(self.fastq_dir, f"{prefix}_{self.in1_suffix}")
        decompress_outfile = os.path.join(self.save_dir, f"{prefix}_{self.out1_suffix}")
        self.check_infile()
        self.check_savedir()
        cmd = (
            f"{self.BBMAP_PROG}"
            f" {self.params}"
            f" {self.infile}"
        )
        super().add_stage("Decompress R1 fastq.gz to fastq", cmd)
        super().add_stage_output_to_file("Write R1 FASTQ file", 0, decompress_outfile, decompress_outfile)

        self.infile = os.path.join(self.fastq_dir, f"{prefix}_{self.in2_suffix}")
        decompress_outfile = os.path.join(self.save_dir, f"{prefix}_{self.out2_suffix}")
        self.check_infile()
        cmd = (
            f"{self.BBMAP_PROG}"
            f" {self.params}"
            f" {self.infile}"
        )
        super().add_stage("Decompress R2 fastq.gz to fastq", cmd)
        super().add_stage_output_to_file("Write R2 FASTQ file", 2, decompress_outfile, decompress_outfile)

    def run(self):
        super().run()
        return all(self.output)


def decompress_demo(config, prefix, fastq_dir="", save_dir=""):
    stage = DecompressStage(config, fastq_dir, save_dir)
    stage.setup(prefix)
    is_complete = stage.run()
    return is_complete