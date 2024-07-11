from edna_processor.analysis_manager import AnalysisManager
from Bio import AlignIO
import numpy as np
import os
import pandas as pd
import tempfile
from edna_processor.utils.utils_hdbscan import fit_hdbscan
from edna_processor.utils.utils_sequence import write_fasta, align_fasta, get_uniq_seq_freq

class HapNetNexusGenerator(AnalysisManager):
 
    def __init__(self, load_path=None):
        super().__init__(load_path)

    def load_points_labels(self,
            index_path: str,
            species_name: str,
            label_type: str
        ) -> list[str]:
        """
        load points labels from index file.
        """
        index = pd.read_csv(index_path, sep='\t')
        subindex = index[index["unit"] == species_name]
        if label_type == 'hdbscan':
            points = subindex[["umap1", "umap2"]].to_numpy()
            labels, _, _, _ = fit_hdbscan(points=points, min_samples=10, min_cluster_size=5)
        elif label_type == 'site':
            labels = ["taoyuan" if "taoyuan" in i else "keelung" for i in subindex["seq_id"]]
        else:
            raise ValueError("Label type must be 'hdbscan' or 'site'")

        return labels

    def write_nexus_file(self,
            species_name: str,
            save_dir: str,
            labels: list[str],
            sample_id_list: list[str]
        ) -> None:
        """
        Write NEXUS file for a given species.
        """
        os.makedirs(save_dir, exist_ok=True)
        temp_dir = tempfile.TemporaryDirectory()

        unit_fasta_path = os.path.join(temp_dir.name, f"{species_name}.fa")
        uniq_unit_fasta_path = os.path.join(temp_dir.name, f"{species_name}_uniq.fa")
        aln_fasta_path = os.path.join(temp_dir.name, f"{species_name}_uniq.aln")
        nex_path = os.path.join(save_dir, f"{species_name}.nex")

        sample_id_list = self.load_sample_id_list(sample_id_list)

        units2fasta = self.load_units2fasta(
            target_name=species_name,
            target_level='species',
            unit_level='species',
            sample_id_list=sample_id_list
        )

        write_fasta(units2fasta, save_path=unit_fasta_path, dereplicate=False)
        write_fasta(units2fasta, save_path=uniq_unit_fasta_path, dereplicate=True)

        freq_string = get_uniq_seq_freq(
            seq_file=unit_fasta_path,
            uniq_seq_file=uniq_unit_fasta_path,
            seq_labels=labels
        )

        align_fasta(seq_file=uniq_unit_fasta_path, aln_file=aln_fasta_path)

        AlignIO.convert(aln_fasta_path, "fasta", nex_path, "nexus", molecule_type="DNA")
        with open(nex_path, 'a') as file:
            file.write(freq_string)

        temp_dir.cleanup()

        print(f"Saved NEXUS file to: {nex_path}")

    def generata_nexus_file(self,
            index_path: str,
            species_name: str,
            label_type: str,
            save_dir: str = '.',
            sample_id_list: list[str] = []
        ) -> None:
        """
        Generate NEXUS file for a given species.
        The NEXUS file contains the aligned sequences, their corresponding labels, and the frequency of each label for each unique sequence.
        'hdbscan' label type uses HDBSCAN clustering results for labeling and calculates the frequency of each cluster label.
        'site' label type uses site information (e.g., 'taoyuan' or 'keelung') for labeling and calculates the frequency of each site label.

        :param index_path: Path to the index file containing the unit information.
        :param species_name: Name of the species for which the NEXUS file will be generated.
        :param label_type: Type of labels to use. Either 'hdbscan' or 'site'.
        :param save_dir: Directory where the NEXUS file will be saved. Default is the current directory.
        :param sample_id_list: List of sample IDs to include in the NEXUS file. The list should be same as that specified by the index file. If empty, all samples will be included.
        """

        labels = self.load_points_labels(
            index_path=index_path,
            species_name=species_name,
            label_type=label_type
        )

        self.write_nexus_file(
            species_name=species_name,
            save_dir=save_dir,
            labels=labels,
            sample_id_list=sample_id_list
        )