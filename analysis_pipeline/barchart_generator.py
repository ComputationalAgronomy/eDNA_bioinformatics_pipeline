from analysis_manager import AnalysisManager
import os
import pandas as pd
from utils_barchart import normalize_abundance, list_union, create_barchart_fig

class BarchartGenerator(AnalysisManager):

    def __init__(self, load_path=None):
        super().__init__(load_path)
        self.samples_abundance = {}

    def get_hap_size(self, sample_id, hap):
        return int(self.sample_data[sample_id].hap_size[hap])

    def get_levelname2abundance_dict(self, sample_id, level):
        abundance = {}
        for hap, level_dict in self.sample_data[sample_id].hap2level.items():
            level_name = level_dict[level]
            if level_name not in abundance:
                abundance[level_name] = 0
            size = self.get_hap_size(sample_id, hap)
            abundance[level_name] += size
        return abundance # e.g. {'SpA': 3, 'SpB': 4, 'SpC': 5}

    def plot_barchart(self, level, save_html_dir=None, save_html_name=None, sample_id_list=None):
        print(f"> Plotting barchart for {level}...")
        sample_id_list = self.load_sample_id_list(sample_id_list)

        for sample_id in sample_id_list:
            abundance = self.get_levelname2abundance_dict(sample_id, level)
            self.samples_abundance[sample_id] = normalize_abundance(abundance)

        all_level_name = [list(self.samples_abundance[sample_id].keys()) for sample_id in sample_id_list]
        uniq_level_name = list_union(all_level_name)

        for sample_id in sample_id_list:
            self.samples_abundance[sample_id] = [self.samples_abundance[sample_id].get(level_name, 0) for level_name in uniq_level_name]

        plotdata = pd.DataFrame(self.samples_abundance, index=uniq_level_name)
        fig = create_barchart_fig(plotdata.transpose())
        fig.show()
        print("> Barchart generated.")

        if save_html_dir is not None:
            save_name = save_name or f"{level}_bar_chart"
            bar_chart_path = os.path.join(save_html_dir, f'{save_html_name}.html')
            fig.write_html(bar_chart_path)
            print(f"> Barchart saved to:  {bar_chart_path}")

        self.analysis_type = "barchart"
        self.results_dir = save_html_dir
        self.parameters.update(
            {
                "level": level,
            }
        )