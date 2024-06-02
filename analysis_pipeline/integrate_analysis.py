from integrate_samples import IntegrateSamples
from run_umap import run_umap
from analysis_function import *
import os
import pandas as pd

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
                if target_name not in rank_dict[target_rank]:
                    continue
                unit_name = rank_dict[units_rank]
                title = f"{unit_name}-{sample_id}_{hap}"
                seq = self.sample_data[sample_id].hap_seq[hap]

                if unit_name not in units2fasta:
                    units2fasta[unit_name] = ""
                units2fasta[unit_name] += f'>{title}\n{seq}\n'
        return units2fasta # e.g. {unit_name: ">SpA_Sample1_Zotu1\nACGT\n>SpA_Sample1_Zotu2\nACGT\n"}

    def mltree_target(self, target_name, target_rank, units_rank="species", dereplicate=False, model=None, bootstrap=None, threads=None, save_dir='.', sample_id_list=None):
        print(f"> Plotting MLTree for {target_name}...")

        ckp = check_mltree_overwrite(save_dir=save_dir, prefix=target_name)
        if ckp == "stop":
            print("> Stopping analysis.")
            return

        create_dir('./tmp/')

        sample_id_list = self.get_sample_id_list(sample_id_list)
        units2fasta = self.get_units2fasta_dict(target_name, target_rank, units_rank, sample_id_list)
        write_fasta(units2fasta, seq_file='./tmp/mltree.fa', dereplicate=dereplicate)
        align_fasta(seq_file='./tmp/mltree.fa', aln_file='./tmp/mltree.aln')
        iqtree2_command(seq_path='./tmp/mltree.aln', save_dir=save_dir, prefix=target_name, model=model, bootstrap=bootstrap, threads=threads, checkpoint=ckp)

        remove_dir('./tmp/')

    def umap_target(self, target_name, target_rank, units_rank="species", neighbors = 15, min_dist = 0.1, save_dir='.', sample_id_list=None):
        print(f"> Plotting UMAP for {target_name}...")

        create_dir('./tmp/')

        sample_id_list = self.get_sample_id_list(sample_id_list)
        units2fasta = self.get_units2fasta_dict(target_name, target_rank, units_rank, sample_id_list)
        write_fasta(units2fasta, seq_file='./tmp/umap.fa', dereplicate=False)
        run_umap(seq_file='./tmp/umap.fa', umap_dir=save_dir, target_name=target_name, neighbors=neighbors, min_dist=min_dist)

        remove_dir('./tmp/')

if __name__ == '__main__':
    read_dir = '.\\data\\all_site'
    load_path = os.path.join(read_dir, 'all_site_sample_data.pkl')
    # a = IntegrateAnalysis()
    # a.import_samples(read_dir=read_dir)
    # a.save_sample_data(save_dir=read_dir, save_name='all_site_sample_data')
    a = IntegrateAnalysis(load_path=load_path)
    # a.umap_target(target_name='Chordata', target_rank='phylum', units_rank='family', neighbors=15, min_dist=0.5, save_dir='.\\test_umap',)
    print(a.sample_id_list)
    # a.barchart_relative_abundance(rank='species')
    # a.umap_target(target_name='Mugil_cephalus', target_rank='species', units_rank='species', neighbors=15, min_dist=0.5, save_dir='.\\test_umap', )
    a.mltree_target(target_name='Mugilidae', target_rank='family', dereplicate=True, model=None, bootstrap=999, threads=None, save_dir='.\\test_mltree')
