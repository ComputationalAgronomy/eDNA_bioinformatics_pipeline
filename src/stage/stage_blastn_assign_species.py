import os
from stage.stage_builder import StageBuilder


class BlastStage(StageBuilder):
    def __init__(self, config, heading="blast", denoise_dir="", save_dir="",
                 db_path: str = "",
                 maxhitnum: int = 1,
                 evalue: float=0.00001,
                 qcov_hsp_perc: int=90,
                 perc_identity: int=90,
                 outfmt: str="10",
                 specifiers: str="qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
        ):
        super().__init__(heading=heading, config=config)
        self.BLAST_PROG = "blastn"
        self.denoise_dir = denoise_dir
        self.save_dir = save_dir
        self.parse_params(db_path, maxhitnum, evalue, qcov_hsp_perc, perc_identity, outfmt, specifiers)

    def parse_params(self, db_path, maxhitnum, evalue, qcov_hsp_perc, perc_identity, outfmt, specifiers):
        self.params = (
            f"-db {db_path}"
            f"-max_target_seqs {maxhitnum}"
            f"-evalue {evalue}"
            f"-qcov_hsp_perc {qcov_hsp_perc}"
            f"-perc_identity {perc_identity}" 
            f'-outfmt "{outfmt} {specifiers}"'
        )

    def setup(self, prefix):
        infile = os.path.join(self.denoise_dir, f"{prefix}_denoise.fasta")
        blast_outfile = os.path.join(self.save_dir, f"{prefix}_blast.csv")
        cmd = (
            f"{self.BLAST_PROG} -query {infile}"
            f"{self.params}"
            f"-num_threads {self.config.n_cpu}"
            f"-out {blast_outfile}"
        )
        super().add_stage("usearch_denoise", cmd)

    def run(self):
        super().run()
        return all(self.output)


def fq_to_fa_demo(config, prefix, denoise_dir="", save_dir=""):
    stage = BlastStage(config, denoise_dir=denoise_dir, save_dir=save_dir)
    stage.setup(prefix)
    is_complete = stage.run()
    return is_complete