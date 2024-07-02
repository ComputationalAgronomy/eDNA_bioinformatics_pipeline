from integrate_samples import IntegrateSamples
from run_umap import write_umap_file, plot_umap
from run_hdbscan import run_hdbscan
from analysis_function import *
import os
import pandas as pd
import numpy as np
from Bio import AlignIO

class IntegrateAnalysis(IntegrateSamples):
    def __init__(self, load_path=None):
        super().__init__(load_path)

    def get_sample_id_list(self, sample_id_list=[]):
        if sample_id_list == []:
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

    def mltree_target(self,
                      target_list, 
                      target_level, 
                      units_level="species", 
                      n_unit_threshold=1, 
                      dereplicate_sequence=True, 
                      model=None, 
                      bootstrap=None, 
                      threads=None, 
                      save_dir='.', 
                      prefix="ml_tree", 
                      sample_id_list=[]):

        print(f"> Plotting MLTree for {" ".join(target_list)}...")

        ckp = check_mltree_overwrite(save_dir=save_dir, prefix=prefix)
        if ckp == "stop":
            print("> Stopping analysis.")
            return

        create_dir('./tmp/')

        sample_id_list = self.get_sample_id_list(sample_id_list)
        fasta_dict = {}
        for target_name in target_list:
            units2fasta = self.get_units2fasta_dict(target_name, target_level, units_level, sample_id_list)
            for unit, fasta in units2fasta.copy().items():
                seq_num = fasta.count('>')
                if seq_num < n_unit_threshold:
                    del units2fasta[unit]
            fasta_dict.update(units2fasta)
        write_fasta(fasta_dict, seq_file='./tmp/mltree.fa', dereplicate=dereplicate_sequence)
        align_fasta(seq_file='./tmp/mltree.fa', aln_file='./tmp/mltree.aln')
        iqtree2_command(seq_path='./tmp/mltree.aln', save_dir=save_dir, prefix=prefix, model=model, bootstrap=bootstrap, threads=threads, checkpoint=ckp)

        remove_dir('./tmp/')

    def umap_target(self,
                    target_list,
                    target_level,
                    units_level="species",
                    index_path=None,
                    neighbors = 15,
                    min_dist = 0.1,
                    dereplicate_sequence=False,
                    parent_dir='.',
                    save_dir_name='umap_result',
                    plot_all=True,
                    plot_target=False,
                    plot_unit=False,
                    plot_target_list=[],
                    cmap="rainbow",
                    sample_id_list=[]):

        print(f"> Plotting UMAP for {" ".join(target_list)}...")

        create_dir('./tmp/')

        if index_path == None:
            sample_id_list = self.get_sample_id_list(sample_id_list)
            fasta_dict = {}
            target2units_dict = {}
            for target_name in target_list:
                units2fasta = self.get_units2fasta_dict(target_name, target_level, units_level, sample_id_list)
                fasta_dict.update(units2fasta)
                target2units_dict[target_name] = list(units2fasta.keys())
            write_fasta(fasta_dict, seq_file='./tmp/umap.fa', dereplicate=dereplicate_sequence)
            save_dir = os.path.join(parent_dir, save_dir_name)
            index = write_umap_file(seq_file='./tmp/umap.fa', save_dir=save_dir, target2units=target2units_dict, neighbors=neighbors, min_dist=min_dist)
        else:
            save_dir = os.path.dirname(index_path)
            index = pd.read_csv(index_path, sep='\t')
        
        index = remove_row_by_unit_occurance(index, n=neighbors)
        if plot_all:
            png_path = os.path.join(save_dir, "all_targets.png")
            plot_umap(index=index, png_path=png_path, cmap=cmap)
        if plot_target:
            targets = index["target"]
            unique_target = np.unique(targets)
            for target in unique_target:
                png_path = os.path.join(save_dir, "family_" + target + ".png")
                subindex = index[index["target"] == target]
                plot_umap(index=subindex, png_path=png_path, cmap=cmap)
        if plot_unit:
            units = index["unit"]
            unique_unit = np.unique(units)
            print(unique_target)
            for unit in unique_unit:
                png_path = os.path.join(save_dir, "species_" + unit + ".png")
                subindex = index[index["unit"] == unit]
                plot_umap(index=subindex, png_path=png_path, cmap=cmap)
        if plot_target_list != [] and type(plot_target_list) == list:
            png_path = os.path.join(save_dir, "family_set.png")
            subindex = index[index["target"].isin(plot_target_list)]
            plot_umap(index=subindex, png_path=png_path, cmap=cmap)

        remove_dir('./tmp/')

    def cluster_umap(self,
                     index_file, 
                     n_unit_threshold=5, 
                     min_samples=5, 
                     min_cluster_size=5, 
                     plot_all=True, 
                     plot_target=False, 
                     plot_unit=False,
                     cmap="rainbow"):

        create_dir("./tmp")        
        cluster_report = []
        save_dir = os.path.dirname(index_file)
        index = pd.read_csv(index_file, sep='\t')
        index = remove_row_by_unit_occurance(index, n_unit_threshold)

        if plot_all:
            png_path = os.path.join(save_dir, "all_targets_clustered.png")
            numb_zotu = len(index)
            numb_unit, numb_clus, clus_perc = run_hdbscan(index, min_samples, min_cluster_size, png_path, cmap)
            cluster_report.append(["all_targets", numb_zotu, numb_unit, numb_clus, clus_perc])
        if plot_target:
            targets = index["target"]
            unique_target = np.unique(targets)
            for target in unique_target:
                png_path = os.path.join(save_dir, "family_" + target + "_clustered.png")
                subindex = index[index["target"] == target]
                units2fasta = self.get_units2fasta_dict(target, "family", "species", self.get_sample_id_list())
                numb_hapl = write_fasta(units2fasta, seq_file='./tmp/cluster.fa', dereplicate=True)
                numb_unit, numb_clus, clus_perc = run_hdbscan(subindex, min_samples, min_cluster_size, png_path, cmap)
                cluster_report.append([target, numb_hapl, numb_unit, numb_clus, clus_perc])
        if plot_unit:
            units = index["unit"]
            unique_unit = np.unique(units)
            for unit in unique_unit:
                png_path = os.path.join(save_dir, "species_" + unit + "_clustered.png")
                subindex = index[index["unit"] == unit]
                units2fasta = self.get_units2fasta_dict(unit, "species", "species", self.get_sample_id_list())
                numb_hapl = write_fasta(units2fasta, seq_file='./tmp/cluster.fa', dereplicate=True)
                numb_unit, numb_clus, clus_perc = run_hdbscan(subindex, min_samples, min_cluster_size, png_path, cmap)
                cluster_report.append([unit, numb_hapl, numb_unit, numb_clus, clus_perc])

        cluster_report_path = os.path.join(save_dir, "cluster_report.tsv")
        cluster_report = pd.DataFrame(cluster_report, columns=["target_name", "number_of_haploids", "number_of_units", "number_of_clusters", "clustered_percentage"])
        cluster_report.to_csv(cluster_report_path, sep='\t', index=False)
        print(f'Saved cluster report to: {cluster_report_path}')
        remove_dir("./tmp")

    def generata_nexus_file(self, index_file, unit_name, save_dir=None):
        create_dir("./tmp")
        if save_dir == None:
            save_dir = os.path.dirname(index_file)
        index = pd.read_csv(index_file, sep='\t')
        subindex = index[index["unit"] == unit_name]
        points = subindex[["umap1", "umap2"]].to_numpy()

        ### HDBSCAN labels
        # labels, _, _, _ = fit_hdbscan(points=points, min_samples=10, min_cluster_size=5)
        ### Site labels color code: #keelung: 44, 126, 247 #taoyuan: 255, 126, 65
        labels = ["taoyuan" if "taoyuan" in i else "keelung" for i in subindex["seq_id"]]

        units2fasta = self.get_units2fasta_dict(unit_name, 'species', 'species', self.get_sample_id_list())
        write_fasta(units2fasta, seq_file=f'./tmp/{unit_name}.fa', dereplicate=False)
        write_fasta(units2fasta, seq_file=f'./tmp/{unit_name}_uniq.fa', dereplicate=True)
        freq_string = get_uniq_seq_freq(seq_file=f'./tmp/{unit_name}.fa', uniq_seq_file=f'./tmp/{unit_name}_uniq.fa', seq_labels=labels)
        align_fasta(seq_file=f'./tmp/{unit_name}_uniq.fa', aln_file=f'./tmp/{unit_name}_uniq.aln')
        AlignIO.convert(f'./tmp/{unit_name}_uniq.aln', "fasta", f"{save_dir}/{unit_name}.nex", "nexus", molecule_type="DNA")
        with open(f"{save_dir}/{unit_name}.nex", 'a') as file:
            file.write(freq_string)
        remove_dir("./tmp")

if __name__ == '__main__':
    read_dir = '.\\data\\all_site'
    load_path = os.path.join(read_dir, 'all_site_alpha_2.pkl')
    # a = IntegrateAnalysis()
    # a.import_samples(read_dir=read_dir)
    # a.save_sample_data(save_dir=read_dir, save_name='amplicon_larger8')
    a = IntegrateAnalysis(load_path=load_path)
    # sample_id_list = a.get_sample_id_list()
    # all_family_name=[]
    # for sample_id in sample_id_list:
    #     for level in a.sample_data[sample_id].hap2level.values():
    #         all_family_name.append(level['family'])
    # uniq_family_name = np.unique(all_family_name)
    # five_families = ["Blenniidae", "Labridae", "Mugilidae", "Pomacentridae", "Tripterygiidae"]
    # a.umap_target(target_list=five_families, target_level='family', dereplicate_sequence=True, units_level='species', neighbors=15, min_dist=0.5, parent_dir='.\\result\\all_site\\permanova', save_dir_name='test', plot_all=False, plot_target=False, plot_unit=False)
    # a.cluster_umap(index_file='.\\..\\..\\result\\all_site\\umap\\five_family_amp\\index.tsv', n_unit_threshold=15, min_samples=30, min_cluster_size=20, plot_all=False, plot_target=True, plot_unit=False)
    # a.barchart_relative_abundance(level='species')
    # a.mltree_target(target_list=five_families, target_level='family', units_level='species', n_unit_threshold=0, dereplicate_sequence=True, bootstrap=1000, threads=None, save_dir='.\\..\\..\\result\\all_site\\mltree\\five_family')
    # a.generata_nexus_file(index_file='.\\..\\..\\result\\all_site\\umap\\five_family_2\\index.tsv', unit_name='Abudefduf_vaigiensis', save_dir=".\\..\\..\\result\\all_site\\hap_net")

    # five_color = {"Blenniidae": "#8000ff", "Labridae": "#00b5eb", "Mugilidae": "#80ffb4", "Pomacentridae": "#ffb360", "Tripterygiidae": "#ff0000"}
    # index = pd.read_csv(".\\result\\all_site\\umap\\five_family_2\\index.tsv", sep='\t')
    # # index = remove_row_by_unit_occurance(index, n=15)
    # target_unit = np.unique(np.array(list(zip(index["target"], index["unit"]))), axis=0)
    # annot_string = ""
    # for i in target_unit:
    #     annot_string += f"CONTAINS=={i[1]} {five_color[i[0]]}\n"
    # print(annot_string)

    index = pd.read_csv(".\\result\\all_site\\umap\\five_family_2\\index.tsv", sep='\t')
    # index = remove_row_by_unit_occurance(index, n=15)
    subindex = index[index["target"] == "Pomacentridae"]
    labels = np.unique(subindex["unit"])
    print(" ".join(labels))
    color = get_color_hex(len(labels))
    annot_string = ""
    for i, label in enumerate(labels):
        annot_string += f"CONTAINS=={label} {color[i]}\n"
        print(i)
    print(annot_string)

    ## itol annot for intraspecific
    # create_dir('./tmp')
    # spc_name = "Enneapterygius_etheostomus"
    # index = pd.read_csv(".\\..\\..\\result\\all_site\\umap\\five_family_2\\index.tsv", sep='\t')
    # subindex = index[index["unit"] == spc_name]
    # points = subindex[["umap1", "umap2"]].to_numpy()
    # labels, _, _, _ = fit_hdbscan(points=points, min_samples=10, min_cluster_size=5)
    # color_code = get_color_hex(len(np.unique(labels)))
    # units2fasta = a.get_units2fasta_dict(spc_name, "species", "species", a.get_sample_id_list())
    # seq_list = re.findall(r"([ATCG]+)\n", units2fasta[spc_name])
    # write_fasta(units2fasta, seq_file='./tmp/seq.fa', dereplicate=True)
    # uniq_records = SeqIO.parse('./tmp/seq.fa', 'fasta')
    # annot_string = ""
    # for rec in uniq_records:
    #     name = rec.id
    #     uniq_seq = str(rec.seq)
    #     for i, seq in enumerate(seq_list):
    #         if uniq_seq==seq:
    #             color = color_code[int(labels[i])]
    #             annot_string += f"CONTAINS=={name} {color}\n"
    #             break
    # print(annot_string)
    # remove_dir('./tmp')
