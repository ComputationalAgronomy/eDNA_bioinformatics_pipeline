from integrate_samples import IntegrateSamples
import os
import pandas as pd
import plotly.express as px

def normalize_abundance(abundance_dict):
    total_size = sum(abundance_dict.values())
    norm_abundance = {key: value/total_size * 100 for key, value in abundance_dict.items()}
    return norm_abundance

def rank_name_union(rank_name_list):
    uniq_rank_name_list = list(set().union(*rank_name_list))
    uniq_rank_name_list.sort()
    return uniq_rank_name_list

def create_barchart_fig(data):
    fig = px.bar(data, barmode='stack', labels={'value': 'Percentage (%)'},color_discrete_sequence=px.colors.qualitative.Pastel)
    fig.update_xaxes(tickmode='linear')
    fig.update_layout(xaxis_title="Sample ID", yaxis_title="Percentage (%)", legend=dict(x=1.05, y=1, traceorder='normal', orientation='h'))
    return fig

class IntegrateAnalysis(IntegrateSamples):
    def __init__(self, load_path=None):
        super().__init__(load_path)

    def barchart_relative_abundance(self, rank, sample_id_list=None, save_dir=None, save_name=None):

        print(f"> Plotting barplot for {rank}...")

        if sample_id_list == None:
            print("> No sample ID list specified. Using all sample IDs.")
            sample_id_list = self.sample_id_list
        else:
            print(f"> Specified sample IDs:  {sample_id_list}")

        all_samples_size = {}
        for sample_id in sample_id_list:
            sample_size = {}
            for hap, rank_list in self.sample_data[sample_id].hap2rank.items():

                rank_name = rank_list[rank]
                size = int(self.sample_data[sample_id].hap_size[hap])

                if rank_name not in sample_size:
                    sample_size[rank_name] = 0
                sample_size[rank_name] += size
            
            all_samples_size[sample_id] = normalize_abundance(sample_size)

        all_rank_name = [list(all_samples_size[sample_id].keys()) for sample_id in sample_id_list]
        uniq_rank_name = rank_name_union(all_rank_name)
    
        for sample_id in sample_id_list:
            all_samples_size[sample_id] = [all_samples_size[sample_id].get(rank_name, 0) for rank_name in uniq_rank_name]
 
        plotdata = pd.DataFrame(all_samples_size, index=uniq_rank_name)
        fig = create_barchart_fig(plotdata.transpose())
        fig.show()
        print("> Barplot generated.")

        if save_dir != None:
            if save_name == None:
                save_name = f"{rank}_bar_chart"
            bar_chart_path = os.path.join(save_dir, f'{save_name}.html')
            fig.write_html(bar_chart_path)
            print(f"> Barplot saved to:  {bar_chart_path}")