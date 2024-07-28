import csv
import numpy as np
import os
import pandas as pd

from fastq_processor.step_build.stage_builder import StageBuilder


class AssignTaxaStage(StageBuilder):

    @staticmethod
    def parse_genus2otherlv(lineage_path: str) -> dict[str, list[str]]:
        genus2otherlv = {}
        with open(lineage_path) as in_handle:
            reader = csv.DictReader(in_handle)
            for row in reader:
                genus_name = row['genus_name']
                otherlv = [
                    row['family_name'],
                    row['order_name'],
                    row['class_name'],
                    row['phylum_name'],
                    row['kingdom_name'],
                ]
                genus2otherlv[genus_name] = otherlv

        return genus2otherlv

    def __init__(self, config, heading="assign_taxa", denoise_dir="", save_dir="",
                 db_path: str = "",
                 lineage_path: str = "",
                 maxhitnum: int = 1,
                 evalue: float = 0.00001,
                 qcov_hsp_perc: int = 90,
                 perc_identity: int = 90,
                 outfmt: str = "10",
                 specifiers: str = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
        ):
        super().__init__(heading=heading, config=config)
        self.BLAST_PROG = "blastn"
        self.in_suffix = "denoise.fasta"
        self.out_suffix = "blast.csv"
        self.denoise_dir = denoise_dir
        self.save_dir = save_dir
        self.blast_outfile = None
        self.lineage_path = lineage_path
        self.parse_params(db_path, maxhitnum, evalue, qcov_hsp_perc, perc_identity, outfmt, specifiers)
        if os.path.isfile(lineage_path):
            self.genus2otherlv = self.parse_genus2otherlv(lineage_path)

    def parse_params(self, db_path, maxhitnum, evalue, qcov_hsp_perc, perc_identity, outfmt, specifiers):
        if db_path == "nt":
            db_path = "nt -remote"
        self.params = (
            f"-db {db_path}"
            f" -max_target_seqs {maxhitnum}"
            f" -evalue {evalue}"
            f" -qcov_hsp_perc {qcov_hsp_perc}"
            f" -perc_identity {perc_identity}" 
            f' -outfmt "{outfmt} {specifiers}"'
        )
        if "remote" not in self.params:
            self.params += f" -num_threads {self.config.n_cpu}"

    def setup(self, prefix):
        infile = os.path.join(self.denoise_dir, f"{prefix}_{self.in_suffix}")
        self.blast_outfile = os.path.join(self.save_dir, f"{prefix}_{self.out_suffix}")
        self.check_path(infile)
        cmd = (
            f"{self.BLAST_PROG} -query {infile}"
            f" {self.params}"
            f" -out {self.blast_outfile}"
        )
        super().add_stage("Taxonomic assignment", cmd)

    def add_taxonomy(self):
        blast_result = pd.read_csv(self.blast_outfile, header=None)
        taxa_matrix = []
        for sseqid in blast_result[1]:
            species = sseqid.split('|')[-1]
            genus = species.split('_')[0]
            if genus in self.genus2otherlv.keys():
                taxonomy_levels = [species, genus] + self.genus2otherlv[genus]
            elif genus.endswith("idae") or genus.endswith("inae"): # only identified to family level
                for _, otherlv in self.genus2otherlv.items():
                    if otherlv[0] == genus:
                        taxonomy_levels = [species, ''] + otherlv
                        break
            elif genus.endswith("iformes") or genus.endswith("oidea"): # only identified to order level
                for _, otherlv in self.genus2otherlv.items():
                    if otherlv[1] == genus:
                        taxonomy_levels = [species, '', ''] + otherlv[1:]
                        break
            else:
                print(f"genus {genus} not found in: {self.lineage_path}")
            taxa_matrix.append(taxonomy_levels)

        taxa_matrix_trans = np.array(taxa_matrix).T
        for i, lv in enumerate(taxa_matrix_trans):
            blast_result.insert(i+2, f"_{i}", lv)

        blast_result.to_csv(self.blast_outfile, index=False, header=None)

        return True

    def run(self):
        super().run()
        if not self.config.dry and hasattr(self, 'genus2otherlv'):
            self.output.append(self.add_taxonomy())
        return all(self.output)


def assign_taxa_demo(config, prefix, denoise_dir="", save_dir="", db_path="", lineage_path="", specifiers="qseqid sscinames sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"):
    stage = AssignTaxaStage(config, denoise_dir=denoise_dir, save_dir=save_dir, db_path=db_path, lineage_path=lineage_path, specifiers=specifiers)
    outfile = stage.setup(prefix)
    is_complete = stage.run()
    return is_complete