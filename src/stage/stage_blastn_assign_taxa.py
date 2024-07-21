import numpy as np
import os
import pandas as pd

from stage.stage_builder import StageBuilder


class AssignTaxaStage(StageBuilder):
    def __init__(self, config, heading="assign_taxa", denoise_dir="", save_dir="",
                 db_path: str = "",
                 maxhitnum: int = 1,
                 evalue: float = 0.00001,
                 qcov_hsp_perc: int = 90,
                 perc_identity: int = 90,
                 outfmt: str = "10",
                 specifiers: str = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
        ):
        super().__init__(heading=heading, config=config)
        self.BLAST_PROG = "blastn"
        self.denoise_dir = denoise_dir
        self.save_dir = save_dir
        self.parse_params(db_path, maxhitnum, evalue, qcov_hsp_perc, perc_identity, outfmt, specifiers)

    def parse_params(self, db_path, maxhitnum, evalue, qcov_hsp_perc, perc_identity, outfmt, specifiers):
        self.params = (
            f"-db {db_path}"
            f" -max_target_seqs {maxhitnum}"
            f" -evalue {evalue}"
            f" -qcov_hsp_perc {qcov_hsp_perc}"
            f" -perc_identity {perc_identity}" 
            f' -outfmt "{outfmt} {specifiers}"'
        )

    def setup(self, prefix):
        self.prefix = prefix
        infile = os.path.join(self.denoise_dir, f"{prefix}_denoise.fasta")
        blast_outfile = os.path.join(self.save_dir, f"{prefix}_blast.csv")
        cmd = (
            f"{self.BLAST_PROG} -query {infile}"
            f" {self.params}"
            f" -num_threads {self.config.n_cpu}"
            f" -out {blast_outfile}"
        )
        super().add_stage("ncbi-blast+_blastn", cmd)

        return blast_outfile

    @staticmethod
    def _add_taxonomy(outfile, genus2taxonomy):
        blast_result = pd.read_csv(outfile, header=None)
        taxa_matrix = []
        for sseqid in blast_result[1]:
            sacc, species = sseqid.split('|')[1], sseqid.split('|')[-1]
            species_firstname = species.split('_')[0]
            if species_firstname in genus2taxonomy.keys():
                genus = species_firstname
            else:
                genus = [genus_level for genus_level, other_levals in genus2taxonomy.items() if species_firstname in other_levals][0]
            taxa_matrix.append([sacc, species, genus, genus2taxonomy[genus][4], genus2taxonomy[genus][3], genus2taxonomy[genus][2], genus2taxonomy[genus][1], genus2taxonomy[genus][0]])
        taxa_matrix = np.array(taxa_matrix)
        sacc_list, species_list, genus_list, family_list, order_list, class_list, phylum_list, kingdom_list = taxa_matrix[:, 0].tolist(), taxa_matrix[:, 1].tolist(), taxa_matrix[:, 2].tolist(), taxa_matrix[:, 3].tolist(), taxa_matrix[:, 4].tolist(), taxa_matrix[:, 5].tolist(), taxa_matrix[:, 6].tolist(), taxa_matrix[:, 7].tolist()
        blast_result[1] = sacc_list
        blast_result.insert(2, '2', species_list)
        blast_result.insert(3, '3', genus_list)
        blast_result.insert(4, '4', family_list)
        blast_result.insert(5, '5', order_list)
        blast_result.insert(6, '6', class_list)
        blast_result.insert(7, '7', phylum_list)
        blast_result.insert(8, '8', kingdom_list)
        blast_result.to_csv(outfile, index=False, header=None)

    def run(self):
        super().run()
        return all(self.output)


def fq_to_fa_demo(config, prefix, denoise_dir="", save_dir=""):
    stage = AssignTaxaStage(config, denoise_dir=denoise_dir, save_dir=save_dir)
    outfile = stage.setup(prefix)
    is_complete = stage.run()
    return is_complete