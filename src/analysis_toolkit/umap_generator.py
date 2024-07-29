import os
import pandas as pd

from analysis_toolkit.data_container import SamplesContainer
from analysis_toolkit.utils.base_logger import logger, get_file_handler
from analysis_toolkit.utils.utils_hdbscan import run_hdbscan_by_category
from analysis_toolkit.utils.utils_umap import run_umap, plot_umap_by_category, filter_index_by_unit_occurrence

class UmapGenerator(SamplesContainer):
    """
    Class for managing UMAP analysis.
    """
    def __init__(self, load_path=None):
        super().__init__(load_path)

    def load_umap_units2fasta(self,
            target_list: list[str],
            target_level: str,
            unit_level: str,
            sample_id_list: list[str]
        ) -> tuple[dict[str, str], dict[str, str]]:
        """
        Updates UMAP units to FASTA mapping for a given set of targets.
        It also creates a dictionary linking unit labels to their corresponding target labels.
        """
        all_targets_units2fasta = {}
        unitlabel2targetlabel = {}
        sample_id_list = self.load_sample_id_list(sample_id_list)

        for target_name in target_list:
            units2fasta = self.load_units2fasta(target_name, target_level, unit_level, sample_id_list)
            all_targets_units2fasta.update(units2fasta)
            unitlabel2targetlabel.update(dict.fromkeys(list(units2fasta.keys()), target_name))

        return all_targets_units2fasta, unitlabel2targetlabel

    def write_umap_index(self,
            target_list: list[str],
            target_level: str, 
            unit_level: str, 
            sample_id_list: list[str], 
            dereplicate_sequence: bool, 
            save_dir: str,
            neighbors: int, 
            min_dist: float, 
            random_state: int, 
            calc_dist: bool
        ) -> pd.DataFrame:
        """
        Write UMAP index to a file and return the index DataFrame.
        """
        units2fasta, unit2target = self.load_umap_units2fasta(
            target_list=target_list,
            target_level=target_level,
            unit_level=unit_level,
            sample_id_list=sample_id_list
        )

        index = run_umap(
            units2fasta=units2fasta,
            dereplicate_sequence=dereplicate_sequence,
            save_dir=save_dir,
            neighbors=neighbors,
            min_dist=min_dist,
            random_state=random_state,
            calc_dist=calc_dist,
            unit2target=unit2target
        )

        return index

    def umap_target(self,
            target_list: list[str],
            target_level: str,
            unit_level: str = "species",
            plot_all: bool = True,
            plot_target: bool = True,
            plot_unit: bool = True,
            save_dir: str = '.',
            n_neighbors: int = 15,
            min_dist: float = 0.1,
            random_state: int = 42,
            calc_dist: bool = True,
            index_path: str = None,
            dereplicate_sequence: bool = False,
            cmap:str = "rainbow",
            show_legend: str = True,
            sample_id_list: list[str] = []
        ) -> None:
        """
        Generate and plot UMAP visualizations for a list of targets based on their sequence data.
        (UMAP parameters reference: https://umap-learn.readthedocs.io/en/latest/parameters.html)

        :param target_list: A list of targets to be plotted (e.g., ["FamilyA", "FamilyB", etc]).
        :param target_level: The taxonomic level of the targets (e.g., family, genus, species).
        :param units_level: The taxonomic level of the units. Default is "species"
        :param plot_all: If True, plot all units in one plot. Default is True
        :param plot_target: If True, plot units for each target in separate plots. Default is True.
        :param plot_unit: If True, plot each unit in separate plots. Default is True
        :param save_dir: Directory where the output files (FASTA, aligned FASTA, index) and plots will be saved. Default is current directory.
        :param n_neighbors: Number of neighbors to consider for UMAP. Default is 15.
        :param min_dist: Minimum distance parameter for UMAP. Default is 0.1.
        :param random_state: Random seed for reproducibility. Default is 42.
        :param calc_dist: If True, calculates a distance matrix for UMAP. Otherwise, transforms sequences into a one-hot encoded matrix. Default is True.
        :param index_path: Path to a precomputed index file. If provided, UMAP embedding calculation will be skipped. Default is None.
        :param dereplicate_sequence: If True, use unique sequences as input data for UMAP. Default is False.
        :param cmap: Colormap for the plots. Default is 'rainbow'.
        :param show_legend: If True, show legend in the plots. Default is True.
        :param sample_id_list: A list of sample IDs to use for UMAP. If not specified (empty), all samples will be used. Default is an empty list.
        """
        os.makedirs(save_dir, exist_ok=True)

        ug_fh = get_file_handler(os.path.join(save_dir, 'umap_generator.log'))
        logger.addHandler(ug_fh)

        logger.info(f"Plotting UMAP for {" ".join(target_list)}...")

        if not index_path:
            index = self.write_umap_index(
                target_list=target_list,
                target_level=target_level,
                unit_level=unit_level,
                sample_id_list=sample_id_list,
                dereplicate_sequence=dereplicate_sequence,
                save_dir=save_dir,
                neighbors=n_neighbors,
                min_dist=min_dist,
                random_state=random_state,
                calc_dist=calc_dist
            )
        else:
            index = pd.read_csv(index_path, sep='\t')

        index = filter_index_by_unit_occurrence(index, n=n_neighbors)

        if plot_all:
            plot_umap_by_category(
                index=index,
                category="all",
                prefix="all_units",
                png_dir=save_dir,
                cmap=cmap,
                show_legend=show_legend
            )

        if plot_target:
            plot_umap_by_category(
                index=index,
                category="target",
                prefix=target_level,
                png_dir=save_dir,
                cmap=cmap,
                show_legend=show_legend
            )

        if plot_unit:
            plot_umap_by_category(
                index=index,
                category="unit",
                prefix=unit_level,
                png_dir=save_dir,
                cmap=cmap,
                show_legend=show_legend
            )

        # record metadata
        self.analysis_type = "UMAP"
        self.results_dir = save_dir
        self.parameters.update(
            {
                "target_list": target_list,
                "target_level": target_level,
                "unit_level": unit_level,
                "plot_all": plot_all,
                "plot_target": plot_target,
                "plot_unit": plot_unit,
                "n_neighbors": n_neighbors,
                "min_dist": min_dist,
                "random_state": random_state,
                "calc_dist": calc_dist,
                "index_path": index_path,
                "dereplicate_sequence": dereplicate_sequence,
                "cmap": cmap,
                "show_legend": show_legend,
            }
        )

    def cluster_umap(self,
            index_file: str,
            n_unit_threshold: int,
            plot_all: bool = False,
            plot_target: bool = False,
            plot_unit: bool = False,
            save_dir: str = '.',
            min_samples: int = 5,
            min_cluster_size: int = 5,
            cluster_selection_epsilon: float = 1.0,
            alpha: float = 1.0,
            cmap: str = "rainbow"
        ) -> None:
        """
        Cluster UMAP embeddings by HDBSCAN and plot the results.
        Clustering could be done at different levels (all, target, unit).
        (HDBSCAN parameters reference: https://hdbscan.readthedocs.io/en/latest/parameter_selection.html)

        :param index_file: Path to the UMAP index file.
        :param n_unit_threshold: Filter out units that occur less than n times in the index DataFrame. The value should be set equal to UMAP n_neighbors.
        :param plot_all: If True, plot all units in one plot. Default is False
        :param plot_target: If True, plot units for each target in separate plots. Default is False.
        :param plot_unit: If True, plot each unit in separate plots. Default is False.
        :param save_dir: Directory where the output cluster report files and plots will be saved. Default is current directory.
        :param min_samples: Provide a measure of how conservative want clustering to be for HDBSCAN. Default is 5.
        :param min_cluster_size: Minimum number of samples required to consider as a cluster for HDBSCAN. Default is 5
        :param cluster_selection_epsilon: The distance threshold that clusters below the given value are not split up any further. Default is 1.0.
        :param alpha: Determines how conservative HDBSCAN will try to cluster points together. Higher values will make HDBSCAN more conservative. Default is 1.0.
        :param cmap: Colormap for the plots.
        """
        os.makedirs(save_dir, exist_ok=True)

        uc_fh = get_file_handler(os.path.join(save_dir, 'umap_cluster.log'))
        logger.addHandler(uc_fh)

        index = pd.read_csv(index_file, sep='\t')
        index = filter_index_by_unit_occurrence(index, n_unit_threshold)

        if plot_all:
            run_hdbscan_by_category(
                index=index,
                min_samples=min_samples,
                min_cluster_size=min_cluster_size,
                cluster_selection_epsilon=cluster_selection_epsilon,
                alpha=alpha,
                category="all",
                prefix="all_units",
                save_dir=save_dir,
                cmap=cmap
            )

        if plot_target:
            run_hdbscan_by_category(
                index=index,
                min_samples=min_samples,
                min_cluster_size=min_cluster_size,
                cluster_selection_epsilon=cluster_selection_epsilon,
                alpha=alpha,
                category="target",
                prefix=self.parameters.get("target_level", "target"),
                save_dir=save_dir,
                cmap=cmap
            )

        if plot_unit:
            run_hdbscan_by_category(
                index=index,
                min_samples=min_samples,
                min_cluster_size=min_cluster_size,
                cluster_selection_epsilon=cluster_selection_epsilon,
                alpha=alpha,
                category="unit",
                prefix=self.parameters.get("unit_level", "unit"),
                save_dir=save_dir,
                cmap=cmap
            )
        
        self.analysis_type = "UMAP"
        self.results_dir = save_dir
        self.parameters.update(
            {
            "index_file": index_file,
            "n_unit_threshold": n_unit_threshold,
            "plot_all": plot_all,
            "plot_target": plot_target,
            "plot_unit": plot_unit,
            "min_samples": min_samples,
            "min_cluster_size": min_cluster_size,
            "cluster_selection_epsilon": cluster_selection_epsilon,
            "alpha": alpha,
            "cmap": cmap,
            }
        )