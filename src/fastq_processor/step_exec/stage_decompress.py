import os
import gzip
from typing import override

from analysis_toolkit.runner_build import base_logger
from fastq_processor.step_build import stage_builder


class DecompressStage(stage_builder.StageBuilder):
    def __init__(self, config, heading="stage_gzip_decompress.py", fastq_dir="", save_dir="",
        ):
        super().__init__(heading=heading, config=config)
        self.insuffix_list = ["R1.fastq.gz", "R2.fastq.gz"]
        self.outsuffix_list = ["R1.fastq", "R2.fastq"]
        self.infile_list = []
        self.outfile_list = []
        self.fastq_dir = fastq_dir
        self.save_dir = save_dir

    def gunzip(infile, outfile):
        with gzip.open(infile, 'rb')as in_handler, open(outfile, 'wb') as out_handler:
            for line in in_handler:
                out_handler.write(line)

    def setup(self, prefix):
        self.check_savedir()
        for in_suffix in self.insuffix_list:
            self.infile = os.path.join(self.fastq_dir, f"{prefix}_{in_suffix}")
            self.check_infile()
            self.infile_list.append(self.infile)

        for out_suffix in self.outsuffix_list:
            self.outfile = os.path.join(self.save_dir, f"{prefix}_{out_suffix}")
            self.outfile_list.append(self.outfile)
        self.prog_name = "Decompress FASTQ.GZ to FASTQ"

    @override
    def run(self):
        base_logger.logger.info(f"Running: {self.heading}")
        base_logger.logger.info(f"Program: {self.prog_name}")
        for infile, outfile in zip(self.infile_list, self.outfile_list):
            DecompressStage.gunzip(infile, outfile)
        base_logger.logger.info(f"COMPLETE: {self.prog_name}")
        return True


def decompress_demo(config, prefix, fastq_dir="", save_dir=""):
    stage = DecompressStage(config, fastq_dir, save_dir)
    stage.setup(prefix)
    is_complete = stage.run()
    return is_complete