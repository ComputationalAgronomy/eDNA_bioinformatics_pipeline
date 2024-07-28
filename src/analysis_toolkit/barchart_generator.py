import os
import pandas as pd

from analysis_toolkit.data_container import SamplesContainer
from analysis_toolkit.utils.base_logger import logger, get_file_handler
from analysis_toolkit.utils.utils_barchart import normalize_abundance, list_union, create_barchart_fig

class BarchartGenerator(SamplesContainer):

    def __init__(self, load_path=None):
        super().__init__(load_path)
        self.samples_abundance = {}

    def load_hap_size(self, sample_id: str, hap: str) -> int:
        """
        Get the size of a haplotype for a given sample.

        :param sample_id: The ID of the sample.
        :param hap: The ID of the haplotype.
        :return: The size of the haplotype.
        """
        return int(self.sample_data[sample_id].hap_size[hap])

    def load_levelname2abundance_dict(self, sample_id: str, level: str):
        """
        Get the abundance of a level for a given sample.

        :param sample_id: The ID of the sample.
        :param level: The name of the level.
        :return: A dictionary mapping level names to abundances.
        """
        abundance = {}
        for hap, level_dict in self.sample_data[sample_id].hap2level.items():
            level_name = level_dict[level]
            if level_name not in abundance:
                abundance[level_name] = 0
            size = self.load_hap_size(sample_id, hap)
            abundance[level_name] += size
        return abundance # e.g. {'SpA': 3, 'SpB': 4, 'SpC': 5}

    def plot_barchart(self,
            level: str,
            save_html_dir: str = '.',
            save_html_name: str = None,
            sample_id_list: list[str] = []
        ) -> None:
        """
        Plot a barchart to visualize the abundance of a level across samples.

        :param level: The name of the level to plot (e.g., species, family, etc.).
        :param save_html_dir: The directory to save the HTML file. Default is current directory.
        :param save_html_name: The name of the HTML file. Default is None.
        :param sample_id_list: A list of sample IDs to plot. Default is None (plot all samples).
        """
        bg_fh = get_file_handler(os.path.join(save_html_dir, "barchart_generator.log"))
        logger.addHandler(bg_fh)

        logger.info(f"Plotting barchart for {level}...")
        sample_id_list = self.load_sample_id_list(sample_id_list)

        for sample_id in sample_id_list:
            abundance = self.load_levelname2abundance_dict(sample_id, level)
            self.samples_abundance[sample_id] = normalize_abundance(abundance)

        all_level_name = [list(self.samples_abundance[sample_id].keys()) for sample_id in sample_id_list]
        uniq_level_name = list_union(all_level_name)

        for sample_id in sample_id_list:
            self.samples_abundance[sample_id] = [self.samples_abundance[sample_id].get(level_name, 0) for level_name in uniq_level_name]

        plotdata = pd.DataFrame(self.samples_abundance, index=uniq_level_name)
        fig = create_barchart_fig(plotdata.transpose())
        fig.show()
        logger.info("Barchart generated.")

        if save_html_dir is not None:
            save_html_name = save_html_name or f"{level}_barchart"
            bar_chart_path = os.path.join(save_html_dir, f'{save_html_name}.html')
            fig.write_html(bar_chart_path)
            logger.info(f"Barchart saved to:  {bar_chart_path}")

        self.analysis_type = "barchart"
        self.results_dir = save_html_dir
        self.parameters.update(
            {
                "level": level,
            }
        )