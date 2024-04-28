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

def list_union(input_list): # e.g. [[1,2,3], [2,3,4], [3,4,5]] -> [1,2,3,4,5]
    uniq_list = list(set().union(*input_list))
    uniq_list.sort()
    return uniq_list

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

def write_umap_fasta(units2seq_dict):
    fasta_path_string = ""
    for unit_name, seq in units2seq_dict.items():
        with open(f'./tmp/{unit_name}.fa', 'w') as file:
            file.write(seq)
        fasta_path_string += f'./tmp/{unit_name}.fa '
    print(f"\n> Written fasta files to: {fasta_path_string}")
    return fasta_path_string # e.g. "./tmp/SpA.fa ./tmp/SpB.fa ./tmp/SpC.fa "

def umap_command(fasta_path_string, save_dir, prefix, neighbors, min_dist):
    save_path = os.path.join(save_dir, prefix)
    cmd = f'usum {fasta_path_string} --neighbors {neighbors} --umap-min-dist {min_dist} --maxdist 1.0 --termdist 1.0 --output {save_path} -f --seed 1'
    subprocess.run(cmd, shell=True)

class IntegrateAnalysis(IntegrateSamples):
    def __init__(self, load_path=None):
        super().__init__(load_path)

    def get_sample_id_list(self, sample_id_list):
        if sample_id_list == None:
            sample_id_list = self.sample_id_list
            print(f"> No sample ID list specified. Using all {len(sample_id_list)} samples.")
        else:
            print(f"> Specified {len(sample_id_list)} samples.")
        return sample_id_list

    def get_hap_size(self, sample_id, hap):
        return int(self.sample_data[sample_id].hap_size[hap])

    def get_rankname2abundance_dict(self, sample_id, rank):
        abundance = {}
        for hap, rank_dict in self.sample_data[sample_id].hap2rank.items():
            rank_name = rank_dict[rank]
            if rank_name not in abundance:
                abundance[rank_name] = 0
            size = self.get_hap_size(sample_id, hap)
            abundance[rank_name] += size
        return abundance # e.g. {'SpA': 3, 'SpB': 4, 'SpC': 5}

    def barchart_relative_abundance(self, rank, save_html_dir=None, save_html_name=None, sample_id_list=None):
        print(f"> Plotting barchart for {rank}...")
        sample_id_list = self.get_sample_id_list(sample_id_list)

        samples_abundance = {}
        for sample_id in sample_id_list:
            abundance = self.get_rankname2abundance_dict(sample_id, rank)
            samples_abundance[sample_id] = normalize_abundance(abundance)

        all_rank_name = [list(samples_abundance[sample_id].keys()) for sample_id in sample_id_list]
        uniq_rank_name = list_union(all_rank_name)
    
        for sample_id in sample_id_list:
            samples_abundance[sample_id] = [samples_abundance[sample_id].get(rank_name, 0) for rank_name in uniq_rank_name]
 
        plotdata = pd.DataFrame(samples_abundance, index=uniq_rank_name)
        fig = create_barchart_fig(plotdata.transpose())
        fig.show()
        print("> Barchart generated.")

        if save_html_dir != None:
            if save_name == None:
                save_name = f"{rank}_bar_chart"
            bar_chart_path = os.path.join(save_html_dir, f'{save_html_name}.html')
            fig.write_html(bar_chart_path)
            print(f"> Barchart saved to:  {bar_chart_path}")

    def get_units2fasta_dict(self, target_name, target_rank, units_rank, sample_id_list):
        units2fasta = {}
        for sample_id in sample_id_list:
            for hap, rank_dict in self.sample_data[sample_id].hap2rank.items():
                if rank_dict[target_rank] != target_name:
                    continue
                unit_name = rank_dict[units_rank]
                title = f"{unit_name}_{sample_id}_{hap}"
                seq = self.sample_data[sample_id].hap_seq[hap]

                if unit_name not in units2fasta:
                    units2fasta[unit_name] = ""
                units2fasta[unit_name] += f'>{title}\n{seq}\n'
        return units2fasta # e.g. {unit_name: ">SpA_Sample1_Zotu1\nACGT\n>SpA_Sample1_Zotu2\nACGT\n"}

    def umap_target(self, target_name, target_rank, units_rank="species", neighbors = 15, min_dist = 0.1, save_dir='.', sample_id_list=None):
        print(f"> Plotting UMAP for {target_name}...")

        create_dir('./tmp/')

        sample_id_list = self.get_sample_id_list(sample_id_list)
        units2fasta = self.get_units2fasta_dict(target_name, target_rank, units_rank, sample_id_list)
        fasta_path = write_umap_fasta(units2fasta)
        umap_command(fasta_path_string=fasta_path, save_dir=save_dir, prefix=target_name, neighbors = neighbors, min_dist = min_dist)

        remove_dir('./tmp/')