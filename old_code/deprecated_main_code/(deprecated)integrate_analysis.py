from integrate_samples import SamplesContainer
from run_hdbscan import run_hdbscan
from Bio import AlignIO
from run_hdbscan import run_hdbscan_by_category
from exec_umap import run_umap, plot_umap_by_category
import numpy as np
import os
import pandas as pd
import tempfile
from utils_barchart import create_barchart_fig, list_union, normalize_abundance
from utils_sequence import write_fasta
from utils_umap import filter_index_by_unit_occurrence
"""
TODO(SW): for IntegrateAnalysis class, which looks like the main class you calling from, this is too complicated. When people look at a function, they can't figure out what does it to within a few seconds. Generally, refactor into separate classes or functions will help.

- i.e., in umap class
def umap_target(self, ...):
    setup(self, ...)
    run_umap(self, ...)
    remove_row_by_unit_occurance(self, ...)
    plot_umap(self, ...)
    cleanup(self, ...)
"""

class AnalysisManager(SamplesContainer):
    def __init__(self, load_path=None):
        super().__init__(load_path)
        self.metadata = {}
        # TODO(SW): My gut feeling is that you should have more variables here, i.e. workspace_dir, target_**

    def load_sample_id_list(self, sample_id_list=[]):
        if sample_id_list == []:
            sample_id_list = self.sample_id_list
            print(f"> No sample ID list specified. Using all {len(sample_id_list)} samples.")
        else:
            print(f"> Specified {len(sample_id_list)} samples.")
        self.metadata["sample_id_list"] = sample_id_list
        return sample_id_list

    def get_units2fasta(self, target_name, target_level, units_level, sample_id_list):
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

class BarchartGenerator(AnalysisManager):
    def __init__(self, load_path=None):
        super().__init__(load_path)
        self.samples_abundance = {}

    def _get_hap_size(self, sample_id, hap):
        return int(self.sample_data[sample_id].hap_size[hap])

    def _get_levelname2abundance_dict(self, sample_id, level):
        abundance = {}
        for hap, level_dict in self.sample_data[sample_id].hap2level.items():
            level_name = level_dict[level]
            if level_name not in abundance:
                abundance[level_name] = 0
            size = self._get_hap_size(sample_id, hap)
            abundance[level_name] += size
        return abundance # e.g. {'SpA': 3, 'SpB': 4, 'SpC': 5}

    def generate_barchart(self, level, save_html_dir=None, save_html_name=None, sample_id_list=None):
        print(f"> Plotting barchart for {level}...")
        sample_id_list = self.load_sample_id_list(sample_id_list)

        self.samples_abundance = {}
        for sample_id in sample_id_list:
            abundance = self._get_levelname2abundance_dict(sample_id, level)
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

        sample_id_list = self.load_sample_id_list(sample_id_list)
        fasta_dict = {}
        for target_name in target_list:
            units2fasta = self.get_units2fasta(target_name, target_level, units_level, sample_id_list)
            for unit, fasta in units2fasta.copy().items():
                seq_num = fasta.count('>')
                if seq_num < n_unit_threshold:
                    del units2fasta[unit]
            fasta_dict.update(units2fasta)
        write_fasta(fasta_dict, seq_file='./tmp/mltree.fa', dereplicate=dereplicate_sequence)
        align_fasta(seq_file='./tmp/mltree.fa', aln_file='./tmp/mltree.aln')
        iqtree2_command(seq_path='./tmp/mltree.aln', save_dir=save_dir, prefix=prefix, model=model, bootstrap=bootstrap, threads=threads, checkpoint=ckp)

        remove_dir('./tmp/')

class UmapAnalyser(AnalysisManager):
    """
    Class for managing UMAP analysis.
    """
    def __init__(self, load_path=None):
        super().__init__(load_path)

    def _update_umap_units2fasta(self,
            target_list: list[str],
            target_level: str,
            units_level: str,
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
            units2fasta = self.get_units2fasta(target_name, target_level, units_level, sample_id_list)
            all_targets_units2fasta.update(units2fasta)
            unitlabel2targetlabel.update(dict.fromkeys(list(units2fasta.keys()), target_name))

        return all_targets_units2fasta, unitlabel2targetlabel

    def _create_umap_index(self,
            target_list: list[str],
            target_level: str, 
            units_level: str, 
            sample_id_list: list[str], 
            dereplicate_sequence: bool, 
            save_dir: str,
            neighbors: int, 
            min_dist: float, 
            random_state: int, 
            calc_dist: bool
        ) -> pd.DataFrame:
        """
        Generate UMAP index for a list of targets.
        """
        temp_dir = tempfile.TemporaryDirectory()
        seq_path = os.path.join(temp_dir.name, 'umap.fa')
    
        units2fasta, unit2target = self._update_umap_units2fasta(
            target_list=target_list,
            target_level=target_level,
            units_level=units_level,
            sample_id_list=sample_id_list
        )
        write_fasta(units2fasta_dict=units2fasta, save_path=seq_path, dereplicate=dereplicate_sequence)

        index = run_umap(
            seq_path=seq_path,
            save_dir=save_dir,
            neighbors=neighbors,
            min_dist=min_dist,
            random_state=random_state,
            calc_dist=calc_dist,
            unit2target=unit2target
        )

        temp_dir.cleanup()

        return index

    def umap_target(self,
            target_list: list[str],
            target_level: str,
            units_level: str = "species",
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
        Generate and plot UMAP visualizations for a list of targets.

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
        print(f"> Plotting UMAP for {" ".join(target_list)}...")
        os.makedirs(save_dir, exist_ok=True)

        if not index_path:
            index = self._create_umap_index(
                target_level=target_list,
                target_level=target_level,
                units_level=units_level,
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
            plot_umap_by_category(index=index, category="all", png_dir=save_dir, cmap=cmap, show_legend=show_legend)
        if plot_target:
            plot_umap_by_category(index=index, category="target", prefix=target_level, png_dir=save_dir, cmap=cmap, show_legend=show_legend)
        if plot_unit:
            plot_umap_by_category(index=index, category="unit", prefix=units_level, png_dir=save_dir, cmap=cmap, show_legend=show_legend)

    def cluster_umap(self,
            index_file,
            plot_all=True,
            plot_target=False,
            plot_unit=False,
            save_dir='.',
            min_samples=5,
            min_cluster_size=5,
            n_unit_threshold=5,
            cmap="rainbow"
        ):
        os.makedirs(save_dir, exist_ok=True)

        index = pd.read_csv(index_file, sep='\t')
        index = filter_index_by_unit_occurrence(index, n_unit_threshold)

        if plot_all:
            cluster_report_all = run_hdbscan_by_category(index, min_samples, min_cluster_size, "all", prefix= save_dir, cmap)
        if plot_target:
            targets = index["target"]
            unique_target = np.unique(targets)
            for target in unique_target:
                png_path = os.path.join(save_dir, "family_" + target + "_clustered.png")
                subindex = index[index["target"] == target]
                units2fasta = self.get_units2fasta(target, "family", "species", self.load_sample_id_list())
                numb_hapl = write_fasta(units2fasta, seq_file='./tmp/cluster.fa', dereplicate=True)
                numb_unit, numb_clus, clus_perc = run_hdbscan(subindex, min_samples, min_cluster_size, png_path, cmap)
                cluster_report.append([target, numb_hapl, numb_unit, numb_clus, clus_perc])
        if plot_unit:
            # TODO(SW): These two blocks are very similar. Can you combine/refactor them? See below.
            units = index["unit"]
            unique_unit = np.unique(units)
            for unit in unique_unit:
                png_path = os.path.join(save_dir, "species_" + unit + "_clustered.png")
                subindex = index[index["unit"] == unit]
                units2fasta = self.get_units2fasta(unit, "species", "species", self.load_sample_id_list())
                numb_hapl = write_fasta(units2fasta, seq_file='./tmp/cluster.fa', dereplicate=True)
                numb_unit, numb_clus, clus_perc = run_hdbscan(subindex, min_samples, min_cluster_size, png_path, cmap)
                cluster_report.append([unit, numb_hapl, numb_unit, numb_clus, clus_perc])

        cluster_report_path = os.path.join(save_dir, "cluster_report.tsv")
        cluster_report = pd.DataFrame(cluster_report, columns=["target_name", "number_of_haploids", "number_of_units", "number_of_clusters", "clustered_percentage"])
        cluster_report.to_csv(cluster_report_path, sep='\t', index=False)
        print(f'Saved cluster report to: {cluster_report_path}')
        remove_dir("./tmp")

    """
    TODO(SW): Refactor thi cluster_report to something like this. Missed a few vars here, please fix that.
    def cluster(self, index, index_key, tax_level, min_samples, min_cluster_size, cmap):
        targets = index[index_key]
        unique_target = np.unique(targets)
        for target in unique_target:
            png_path = os.path.join(save_dir, f"{tax_level}_{target}_clustered.png")
            subindex = index[index[index_key] == target]
            units2fasta = self.get_units2fasta_dict(target, tax_level, "species", self.get_sample_id_list())
            numb_hapl = write_fasta(units2fasta, seq_file='./tmp/cluster.fa', dereplicate=True)
            numb_unit, numb_clus, clus_perc = run_hdbscan(subindex, min_samples, min_cluster_size, png_path, cmap)
            cluster_report.append([target, numb_hapl, numb_unit, numb_clus, clus_perc])

    cluster("target", "family")
    cluster("unit", "species")

    """


    def generata_nexus_file(self, index_file, unit_name, save_dir=None):
        create_dir("./tmp")

        save_dir = save_dir or os.path.dirname(index_file)
        index = pd.read_csv(index_file, sep='\t')
        subindex = index[index["unit"] == unit_name]
        points = subindex[["umap1", "umap2"]].to_numpy()

        ### HDBSCAN labels
        # labels, _, _, _ = fit_hdbscan(points=points, min_samples=10, min_cluster_size=5)
        ### Site labels color code: #keelung: 44, 126, 247 #taoyuan: 255, 126, 65
        labels = ["taoyuan" if "taoyuan" in i else "keelung" for i in subindex["seq_id"]]

        units2fasta = self.get_units2fasta(unit_name, 'species', 'species', self.load_sample_id_list())
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
    sample_id_list = a.get_sample_id_list()
    all_family_name=[]
    for sample_id in sample_id_list:
        for level in a.sample_data[sample_id].hap2level.values():
            all_family_name.append(level['family'])
    uniq_family_name = np.unique(all_family_name)
    # five_families = ["Blenniidae", "Labridae", "Mugilidae", "Pomacentridae", "Tripterygiidae"]
    a.umap_target(target_list=uniq_family_name, target_level='family', dereplicate_sequence=False, units_level='species', calc_dist=False, n_neighbors=15, min_dist=0.5, parent_dir='.\\result\\all_site\\umap', save_dir_name='all_family_test_number_ver', plot_all=True, plot_target=True, plot_unit=True)
    # a.cluster_umap(index_file='.\\result\\all_site\\umap\\all_family_test_number_ver\\index.tsv', n_unit_threshold=15, min_samples=10, min_cluster_size=10, plot_all=False, plot_target=False, plot_unit=True)
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

    # index = pd.read_csv(".\\result\\all_site\\umap\\five_family_2\\index.tsv", sep='\t')
    # # index = remove_row_by_unit_occurance(index, n=15)
    # subindex = index[index["target"] == "Pomacentridae"]
    # labels = np.unique(subindex["unit"])
    # print(" ".join(labels))
    # color = get_color_hex(len(labels))
    # annot_string = ""
    # for i, label in enumerate(labels):
    #     annot_string += f"CONTAINS=={label} {color[i]}\n"
    #     print(i)
    # print(annot_string)

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
