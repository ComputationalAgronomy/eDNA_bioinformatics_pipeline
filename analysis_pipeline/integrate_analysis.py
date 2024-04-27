from integrate_samples import IntegrateSamples
import os
import pandas as pd
import plotly.express as px
import subprocess
import shutil

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

def create_dir(dir_path):
    if not os.path.isdir(dir_path):
        print(f"> Creating directory: {dir_path}")
        os.makedirs(dir_path)

def remove_dir(dir_path):
    if os.path.isdir(dir_path):
        print(f"> Removing directory: {dir_path}")
        shutil.rmtree(dir_path)

def umap_command(files_string, save_path, neighbors = 15, min_dist = 0.1):
    cmd = f'usum {files_string} --neighbors {neighbors} --umap-min-dist {min_dist} --maxdist 1.0 --termdist 1.0 --output {save_path} -f --seed 1'
    subprocess.run(cmd, shell=True)

class IntegrateAnalysis(IntegrateSamples):
    def __init__(self, load_path=None):
        super().__init__(load_path)

    def get_sample_abundance(self, sample_id, rank):

        abundance = {}
        for hap, rank_dict in self.sample_data[sample_id].hap2rank.items():
            rank_name = rank_dict[rank]
            if rank_name not in abundance:
                abundance[rank_name] = 0
            size = int(self.sample_data[sample_id].hap_size[hap])
            abundance[rank_name] += size

        return abundance

    def barchart_relative_abundance(self, rank, sample_id_list=None, save_dir=None, save_name=None):

        print(f"> Plotting barchart for {rank}...")

        if sample_id_list == None:
            print("> No sample ID list specified. Using all samples.")
            sample_id_list = self.sample_id_list
        else:
            print(f"> Specified samples:  {sample_id_list}")

        samples_abundance = {}
        for sample_id in sample_id_list:
            abundance = self.get_sample_abundance(sample_id, rank)
            samples_abundance[sample_id] = normalize_abundance(abundance)

        all_rank_name = [list(samples_abundance[sample_id].keys()) for sample_id in sample_id_list]
        uniq_rank_name = rank_name_union(all_rank_name)
    
        for sample_id in sample_id_list:
            samples_abundance[sample_id] = [samples_abundance[sample_id].get(rank_name, 0) for rank_name in uniq_rank_name]
 
        plotdata = pd.DataFrame(samples_abundance, index=uniq_rank_name)
        fig = create_barchart_fig(plotdata.transpose())
        fig.show()
        print("> Barchart generated.")

        if save_dir != None:
            if save_name == None:
                save_name = f"{rank}_bar_chart"
            bar_chart_path = os.path.join(save_dir, f'{save_name}.html')
            fig.write_html(bar_chart_path)
            print(f"> Barchart saved to:  {bar_chart_path}")

    def umap_target(self, target_name, target_rank, units_rank="species", save_dir='.', neighbors = 15, min_dist = 0.1, sample_id_list=None):

        create_dir('./tmp/')

        print(f"> Plotting UMAP for {target_name}...")

        if sample_id_list == None:
            print("> No sample ID list specified. Using all samples.")
            sample_id_list = self.sample_id_list
        else:
            print(f"> Specified samples:  {sample_id_list}")

        units2seq = {}
        for sample_id in sample_id_list:
            for hap, rank_dict in self.sample_data[sample_id].hap2rank.items():
                if rank_dict[target_rank] != target_name:
                    continue
                unit_name = rank_dict[units_rank]
                title = f"{unit_name}_{sample_id}_{hap}"
                seq = self.sample_data[sample_id].hap_seq[hap]

                if unit_name not in units2seq:
                    units2seq[unit_name] = ""
                units2seq[unit_name] += f'>{title}\n{seq}\n'

        file_string = ""
        for unit_name, seq in units2seq.items():
            with open(f'./tmp/{unit_name}.fasta', 'w') as file:
                file.write(seq)
            file_string = file_string + f'./tmp/{unit_name}.fasta '

        save_path = os.path.join(save_dir, target_name)
        umap_command(files_string=file_string, save_path=save_path, neighbors = neighbors, min_dist = min_dist)

        remove_dir('./tmp/')
