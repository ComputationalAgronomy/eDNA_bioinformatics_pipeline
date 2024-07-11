from edna_processor.analysis_manager import AnalysisManager
import os
import tempfile
from edna_processor.utils_sequence import write_fasta, align_fasta
from edna_processor.utils_tree import run_iqtree2

class MLTreeGenerator(AnalysisManager):

    def __init__(self, load_path=None):
        super().__init__(load_path)

    def load_mltree_units2fasta(self,
            target_list: list[str],
            target_level: str,
            unit_level: str,
            n_unit_threshold: int,
            sample_id_list: list[str]
        ) -> dict[str, str]:
        """
        Load units2fasta dictionary for a list of targets.
        """
        fasta_dict = {}
        for target_name in target_list:
            units2fasta = self.load_units2fasta(target_name, target_level, unit_level, sample_id_list)

            for unit, fasta in units2fasta.copy().items():
                seq_num = fasta.count('>')
                if seq_num < n_unit_threshold:
                    del units2fasta[unit]

            fasta_dict.update(units2fasta)

        return fasta_dict

    def write_mltree_fasta(self,
            target_list: list[str],
            target_level: str,
            unit_level: str,
            n_unit_threshold: int,
            dereplicate_sequence: bool,
            save_path: str,
            sample_id_list: list[str]
        ) -> None:
        """
        Write an aligned FASTA file for a list of targets.
        """
        unit2fasta = self.load_mltree_units2fasta(
            target_list = target_list,
            target_level=target_level,
            unit_level=unit_level,
            n_unit_threshold=n_unit_threshold,
            sample_id_list=sample_id_list
        )
        
        temp_dir = tempfile.TemporaryDirectory()

        unit_fasta_path = os.path.join(temp_dir.name, 'mltree.fa')

        write_fasta(unit2fasta, save_path=unit_fasta_path, dereplicate=dereplicate_sequence)
        align_fasta(seq_file=unit_fasta_path, aln_file=save_path)

        temp_dir.cleanup()

    def mltree_target(self,
            target_list: list[str],
            target_level: str,
            unit_level:str = "species",
            save_dir:str = '.',
            prefix: str = "ml_tree",
            model: str = None,
            bootstrap: int = None,
            threads: int = None,
            dereplicate_sequence: bool = True,
            n_unit_threshold:int = 1,
            sample_id_list: list[str] = []
        ) -> None:
        """
        Reconstruct a phylogenetic tree for a list of targets using IQTREE.
        (IQTREE2 command reference: http://www.iqtree.org/doc/Command-Reference)

        :param target_list: A list of targets to be plotted (e.g., ["FamilyA", "FamilyB", etc]).
        :param target_level: The taxonomic level of the targets (e.g., family, genus, species).
        :param units_level: The taxonomic level of the units. Default is "species"
        :param save_dir: Directory to save the output files.
        :param prefix: Prefix for the output file names.
        :param model: Model to specify for tree inference. If not specified, it will use the best-fit model found. Default is None.
        :param bootstrap: Number of bootstrap replicates. Default is None.
        :param threads: Number of threads to use. If not specified, it will automatically determine the best number of cores given the current data and computer. Default is None.
        """

        print(f"> Plotting MLTree for {" ".join(target_list)}...")

        os.makedirs(save_dir, exist_ok=True)
 
        ml_fasta_path = os.path.join(save_dir, f'{prefix}.aln')

        sample_id_list = self.load_sample_id_list(sample_id_list)

        self.write_mltree_fasta(
            target_list=target_list,
            target_level=target_level,
            unit_level=unit_level,
            n_unit_threshold=n_unit_threshold,
            dereplicate_sequence=dereplicate_sequence,
            save_path=ml_fasta_path,
            sample_id_list=sample_id_list
        )

        run_iqtree2(
            seq_path=ml_fasta_path,
            save_dir=save_dir,
            prefix=prefix,
            model=model,
            bootstrap=bootstrap,
            threads=threads
        )

        self.analysis_type = "ml_tree"
        self.results_dir = save_dir
        self.parameters.update(
            {
                "target_list": target_list,
                "target_level": target_level,
                "unit_level": unit_level,
                "prefix": prefix,
                "model": model,
                "bootstrap": bootstrap,
                "threads": threads,
                "dereplicate_sequence": dereplicate_sequence,
                "n_unit_threshold": n_unit_threshold,
                "sample_id_list": sample_id_list
            }
        )