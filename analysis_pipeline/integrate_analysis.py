from integrate_samples import IntegrateSamples
from run_umap import write_umap_file, plot_umap
from run_hdbscan import run_hdbscan
from analysis_function import *
import os
import pandas as pd
import numpy as np

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

    def get_levelname2abundance_dict(self, sample_id, level):
        abundance = {}
        for hap, level_dict in self.sample_data[sample_id].hap2level.items():
            level_name = level_dict[level]
            if level_name not in abundance:
                abundance[level_name] = 0
            size = self.get_hap_size(sample_id, hap)
            abundance[level_name] += size
        return abundance # e.g. {'SpA': 3, 'SpB': 4, 'SpC': 5}

    def barchart_relative_abundance(self, level, save_html_dir=None, save_html_name=None, sample_id_list=None):
        print(f"> Plotting barchart for {level}...")
        sample_id_list = self.get_sample_id_list(sample_id_list)

        samples_abundance = {}
        for sample_id in sample_id_list:
            abundance = self.get_levelname2abundance_dict(sample_id, level)
            samples_abundance[sample_id] = normalize_abundance(abundance)

        all_level_name = [list(samples_abundance[sample_id].keys()) for sample_id in sample_id_list]
        uniq_level_name = list_union(all_level_name)
    
        for sample_id in sample_id_list:
            samples_abundance[sample_id] = [samples_abundance[sample_id].get(level_name, 0) for level_name in uniq_level_name]
 
        plotdata = pd.DataFrame(samples_abundance, index=uniq_level_name)
        fig = create_barchart_fig(plotdata.transpose())
        fig.show()
        print("> Barchart generated.")

        if save_html_dir != None:
            if save_name == None:
                save_name = f"{level}_bar_chart"
            bar_chart_path = os.path.join(save_html_dir, f'{save_html_name}.html')
            fig.write_html(bar_chart_path)
            print(f"> Barchart saved to:  {bar_chart_path}")

    def get_units2fasta_dict(self, target_name, target_level, units_level, sample_id_list):
        units2fasta = {}
        for sample_id in sample_id_list:
            for hap, level_dict in self.sample_data[sample_id].hap2level.items():
                if target_name not in level_dict[target_level]:
                    continue
                unit_name = level_dict[units_level]
                title = f"{unit_name}-{sample_id}_{hap}"
                seq = self.sample_data[sample_id].hap_seq[hap]

                if unit_name not in units2fasta:
                    units2fasta[unit_name] = ""
                units2fasta[unit_name] += f'>{title}\n{seq}\n'
        return units2fasta # e.g. {unit_name: ">SpA_Sample1_Zotu1\nACGT\n>SpA_Sample1_Zotu2\nACGT\n"}

    def mltree_target(self, target_name, target_level, units_level="species", dereplicate=False, model=None, bootstrap=None, threads=None, save_dir='.', sample_id_list=None):
        print(f"> Plotting MLTree for {target_name}...")

        ckp = check_mltree_overwrite(save_dir=save_dir, prefix=target_name)
        if ckp == "stop":
            print("> Stopping analysis.")
            return

        create_dir('./tmp/')

        sample_id_list = self.get_sample_id_list(sample_id_list)
        units2fasta = self.get_units2fasta_dict(target_name, target_level, units_level, sample_id_list)
        write_fasta(units2fasta, seq_file='./tmp/mltree.fa', dereplicate=dereplicate)
        align_fasta(seq_file='./tmp/mltree.fa', aln_file='./tmp/mltree.aln')
        iqtree2_command(seq_path='./tmp/mltree.aln', save_dir=save_dir, prefix=target_name, model=model, bootstrap=bootstrap, threads=threads, checkpoint=ckp)

        remove_dir('./tmp/')

    def umap_target(self,
                    target_list,
                    target_level,
                    units_level="species",
                    neighbors = 15,
                    min_dist = 0.1,
                    parent_dir='.',
                    save_dir_name='umap_result',
                    plot_all=True,
                    plot_sep=False,
                    cmap="rainbow",
                    sample_id_list=None):

        print(f"> Plotting UMAP for {" ".join(target_list)}...")

        create_dir('./tmp/')

        sample_id_list = self.get_sample_id_list(sample_id_list)
        fasta_dict = {}
        target2units_dict = {}
        for target_name in target_list:
            units2fasta = self.get_units2fasta_dict(target_name, target_level, units_level, sample_id_list)
            fasta_dict.update(units2fasta)
            target2units_dict[target_name] = list(units2fasta.keys())
        write_fasta(fasta_dict, seq_file='./tmp/umap.fa', dereplicate=False)
        save_dir = os.path.join(parent_dir, save_dir_name)
        index = write_umap_file(seq_file='./tmp/umap.fa', save_dir=save_dir, target2units=target2units_dict, neighbors=neighbors, min_dist=min_dist)

        if plot_all:
            png_path = os.path.join(save_dir, "umap_all_targets.png")
            plot_umap(index=index, n_unit_threshold=neighbors, png_path=png_path, cmap=cmap)
        if plot_sep:
            targets = index["target"]
            unique_target = np.unique(targets)
            for target in unique_target:
                png_path = os.path.join(save_dir, target + ".png")
                subindex = index[index["target"] == target]
                plot_umap(index=subindex, n_unit_threshold=neighbors, png_path=png_path, cmap=cmap)

        remove_dir('./tmp/')

    @staticmethod
    def cluster_umap(index_file, 
                     n_unit_threshold=5, 
                     min_samples=1, 
                     min_cluster_size=5, 
                     plot_all=True, 
                     plot_sep=False, 
                     cmap="rainbow"):
        
        cluster_report = []
        index = pd.read_csv(index_file, sep='\t')
        save_dir = os.path.dirname(index_file)

        if plot_all:
            png_path = os.path.join(save_dir, "all_targets_clustered.png")
            numb_unit, numb_clus, clus_perc = run_hdbscan(index, n_unit_threshold, min_samples, min_cluster_size, png_path, cmap)
            cluster_report.append(["all_targets", numb_unit, numb_clus, clus_perc])
        if plot_sep:
            targets = index["target"]
            unique_target = np.unique(targets)
            for target in unique_target:
                png_path = os.path.join(save_dir, target + "_clustered.png")
                subindex = index[index["target"] == target]
                numb_unit, numb_clus, clus_perc = run_hdbscan(subindex, n_unit_threshold, min_samples, min_cluster_size, png_path, cmap)
                cluster_report.append([target, numb_unit, numb_clus, clus_perc])

        cluster_report_path = os.path.join(save_dir, "cluster_report.tsv")
        cluster_report = pd.DataFrame(cluster_report, columns=["target_name", "number_of_units", "number_of_clusters", "clustered_percentage"])
        cluster_report.to_csv(cluster_report_path, sep='\t', index=False)
        print(f'Saved cluster report to: {cluster_report_path}')

if __name__ == '__main__':
    read_dir = '.\\..\\..\\data\\all_site'
    load_path = os.path.join(read_dir, 'all_site_alpha_2.pkl')
    # a = IntegrateAnalysis()
    # a.import_samples(read_dir=read_dir)
    # a.save_sample_data(save_dir=read_dir, save_name='all_site_alpha_8')
    a = IntegrateAnalysis(load_path=load_path)
    # a.umap_target(target_name='Chordata', target_level='phylum', units_level='family', neighbors=15, min_dist=0.5, save_dir='.\\test_umap',)
    # print(a.sample_id_list)
    # a.barchart_relative_abundance(level='species')
    # a.umap_target(target_list=["Blenniidae", "Labridae", "Mugilidae", "Pomacentridae", "Tripterygiidae"], target_level='family', units_level='species', neighbors=15, min_dist=0.5, parent_dir='.\\..\\..\\test_umap', plot_all=True, plot_sep=True,)
    a.cluster_umap(index_file='.\\..\\..\\test_umap\\umap_result\\index.tsv', n_unit_threshold=15, min_samples=5, min_cluster_size=5, plot_all=True, plot_sep=True)
    # a.mltree_target(target_name='Mugilidae', target_level='family', dereplicate=True, model=None, bootstrap=999, threads=None, save_dir='.\\test_mltree')
    # a.umap_target(target_name='Mugilidae', target_level='family', units_level='species', neighbors=15, min_dist=0.5, save_dir='.\\test_umap')
