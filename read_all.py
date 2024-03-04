from abc import ABC, abstractmethod
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import subprocess
import re
from pymsaviz import MsaViz
import shutil
import plotly.express as px
from usum import usum

class OtuAnalysis(ABC):

    def __init__(self, read_dir, sample_num):
        self.uniq_seq = {}
        self.uniq_size = {}
        self.otu_seq = {}
        self.otu_size = {}
        self.otu2unique = {}
        self.taxonomy2otu = {'kingdom':{}, 'phylum':{}, 'class':{}, 'order':{}, 'family':{}, 'genus':{}, 'species':{}}
        self.taxonomylevel= {'kingdom':{}, 'phylum':{}, 'class':{}, 'order':{}, 'family':{}, 'genus':{}, 'species':{}}
        uniq_path = f'{read_dir}/4_derep/{sample_num}_derep.fasta'
        otu_path = f'{read_dir}/5_haploid/{sample_num}_zotu.fasta'
        otu_report_path = f'{read_dir}/5_haploid/{sample_num}_zotu_report.txt'
        blast_path = f'{read_dir}/6_blast/{sample_num}_zotu.csv'
        
        self.uniq_seq = self._read_seq(uniq_path)
        self.otu_seq = self._read_seq(otu_path)
        self._read_otu_report(otu_report_path)
        self._get_otu_size()
        self._read_blast(blast_path)

    def _read_seq(self, seq_path, prefix = '>'):
        with open(seq_path,'r') as file:
            seq={}
            haploid_pattern = re.compile(f'{prefix}[a-zA-Z]+\d+')
            for line in file.readlines():
                match = haploid_pattern.match(line)
                if match:
                    haploid = match.group().replace(prefix,'')
                    seq[haploid] = ''
                else:
                    seq[haploid] = seq[haploid] + line.strip()
            return seq

    @abstractmethod
    def _read_otu_report(report_path):
        raise NotImplementedError

    def _get_otu_size(self):
        for otu, uniq_group in self.otu2unique.items():
            size = sum(int(self.uniq_size[uniq]) for uniq in uniq_group)
            self.otu_size[otu] = size

    def _read_blast(self, blast_path):
        with open(blast_path, 'r') as file:
            self.taxonomy2otu['order']['no_order_level'] = []
            self.taxonomylevel['order']['no_order_level'] = []

            for line in file.readlines():
                blast_list = line.split(',')
                haploid, species, genus, family, order, class_, phylum, kingdom = blast_list[0], blast_list[2], blast_list[3], blast_list[4], blast_list[5], blast_list[6], blast_list[7], blast_list[8]
                error_symbol= str.maketrans({':': '_', '/': '_', '\\': '_', '*': '_', '?': '_', '"': '_', '<': '_', '>': '_', '|': '_'})
                species = species.translate(error_symbol)
                family = 'Mugilidae' if family == 'Mugil' else family
            
                if species not in self.taxonomy2otu['species']:
                    self.taxonomy2otu['species'][species] = []
                self.taxonomy2otu['species'][species].append(haploid)

                if genus not in self.taxonomy2otu['genus']:
                    self.taxonomy2otu['genus'][genus] = []
                self.taxonomy2otu['genus'][genus].append(haploid)

                if family not in self.taxonomy2otu['family']:
                    self.taxonomy2otu['family'][family] = []
                self.taxonomy2otu['family'][family].append(haploid)

                if order == '':
                    self.taxonomy2otu['order']['no_order_level'].append(haploid)
                else:
                    if order not in self.taxonomy2otu['order']:
                        self.taxonomy2otu['order'][order] = []
                    self.taxonomy2otu['order'][order].append(haploid)

                if class_ not in self.taxonomy2otu['class']:
                    self.taxonomy2otu['class'][class_] = []
                self.taxonomy2otu['class'][class_].append(haploid)

                if phylum not in self.taxonomy2otu['phylum']:
                    self.taxonomy2otu['phylum'][phylum] = []
                self.taxonomy2otu['phylum'][phylum].append(haploid)

                if kingdom not in self.taxonomy2otu['kingdom']:
                    self.taxonomy2otu['kingdom'][kingdom] = []
                self.taxonomy2otu['kingdom'][kingdom].append(haploid)


                if species not in self.taxonomylevel['species']:
                    self.taxonomylevel['species'][species] = []
                self.taxonomylevel['species'][species].append(haploid)

                if genus not in self.taxonomylevel['genus']:
                    self.taxonomylevel['genus'][genus] = []
                if species not in self.taxonomylevel['genus'][genus]:
                    self.taxonomylevel['genus'][genus].append(species)

                if family not in self.taxonomylevel['family']:
                    self.taxonomylevel['family'][family] = []
                if genus not in self.taxonomylevel['family'][family]:
                    self.taxonomylevel['family'][family].append(genus)

                if order == '' and family not in self.taxonomylevel['order']['no_order_level']:
                    self.taxonomylevel['order']['no_order_level'].append(family)
                else:
                    if order not in self.taxonomylevel['order']:
                        self.taxonomylevel['order'][order] = []
                    if family not in self.taxonomylevel['order'][order]:
                        self.taxonomylevel['order'][order].append(family)

                if class_ not in self.taxonomylevel['class']:
                    self.taxonomylevel['class'][class_] = []
                if order not in self.taxonomylevel['class'][class_]:
                    self.taxonomylevel['class'][class_].append(order)

                if phylum not in self.taxonomylevel['phylum']:
                    self.taxonomylevel['phylum'][phylum] = []
                if class_ not in self.taxonomylevel['phylum'][phylum]:
                    self.taxonomylevel['phylum'][phylum].append(class_)

                if kingdom not in self.taxonomylevel['kingdom']:
                    self.taxonomylevel['kingdom'][kingdom] = []
                if phylum not in self.taxonomylevel['kingdom'][kingdom]:
                    self.taxonomylevel['kingdom'][kingdom].append(phylum)

    # Split the sequences to conform to the fasta format.
    @staticmethod
    def _split_lines(seq):
        lines = [seq[i:i+59] for i in range(0, len(seq), 60)]
        return '\n'.join(lines)
    
    @staticmethod
    def _make_tmp_dir():
        if not os.path.isdir('./tmp/'):
            os.makedirs('./tmp/')

    @staticmethod
    def _make_umap(file_string, save_folder_name, neighbors = 15, umap_min_dist = 0.1):
        cmd = f'usum {file_string} --neighbors {neighbors} --umap-min-dist {umap_min_dist} --maxdist 1.0 --termdist 1.0 --output {save_folder_name} -f'
        subprocess.run(cmd, shell=True)

    @staticmethod
    def _make_align_file(seq_file, aln_file):
        cmd = f'clustalo -i {seq_file} -o {aln_file}.aln --distmat-out={aln_file}.mat --guidetree-out={aln_file}.dnd --full --force'
        subprocess.run(cmd, shell=True)

    @staticmethod
    def _show_alignment(aln_seq_path, save_path, save_png=True):
        mv = MsaViz(aln_seq_path, wrap_length=60, show_count=True, show_grid=False, show_consensus=True)
        mv.plotfig()
        if save_png==True:
            mv.savefig(save_path)

    #plot MSA between unique sequences in one OTU.
    def within_otu_align(self, otu_name, save_dir='.', save_align_png=True):
        self._make_tmp_dir()

        uniq_list = self.otu2unique[otu_name]
        text = ""
        for uniq in uniq_list:
            text = text + f'>{uniq}\n{self.uniq_seq[uniq]}\n'
        with open('./tmp/seq.fa', 'w') as f:
            f.write(text)
        
        self._make_align_file(seq_file='./tmp/seq.fa', aln_file='./tmp/seq.aln')

        if os.path.exists('./tmp/seq.aln'):
            self._show_alignment(aln_seq_path='./tmp/seq.aln', save_path=f'{save_dir}/{otu_name}.png', save_png=save_align_png)
        else:
            print(f'{otu_name} contains 1 sequence, nothing to align')
        
        shutil.rmtree('./tmp/')

    # plot UMAP and alignment between OTUs in one species.
    def analysis_species(self, species_name, save_dir='.', save_align_png=True, make_phylogenetic_tree=False, bootstrap_times=100):
        self._make_tmp_dir()

        otu_list = self.taxonomy2otu['species'][species_name]
        for otu_name in otu_list:
            seq = ''
            for uniq in self.otu2unique[otu_name]:
                subseq = self.uniq_seq[uniq]
                subseq = self._split_lines(subseq)
                seq = seq + f'>{uniq}\n{subseq}\n'
            with open(f'./tmp/{otu_name}.fasta', 'w') as f:
                f.write(seq)
        file_list = ['./tmp/' + otu + '.fasta' for otu in otu_list]
        file_string = ' '.join(file_list)
        self._make_umap(file_string=file_string, save_folder_name=f'{save_dir}/{species_name}')

        seq = ''
        for otu in otu_list:
            subseq = self.otu_seq[otu]
            subseq = self._split_lines(subseq)
            seq = seq + f'>{otu}\n{subseq}\n'
        with open('./tmp/seq.fa', 'w') as f:
            f.write(seq)
        
        self._make_align_file(seq_file='./tmp/seq.fa', aln_file='./tmp/seq.aln')
        
        if os.path.exists('./tmp/seq.aln'):
            self._show_alignment(aln_seq_path='./tmp/seq.aln', save_path=f'{save_dir}/{species_name}/alignment.png', save_png=save_align_png)
        else:
            print(f'{species_name} contains 1 sequence, nothing to align')

        if make_phylogenetic_tree==True:
            cmd = f'iqtree2 -m GTR+F+I+G4 -s ./tmp/seq.aln -b {bootstrap_times} --prefix {species_name}'
            subprocess.run(cmd, shell=True)
            cmd = f'move {species_name}.* {save_dir}/{species_name}'
            subprocess.run(cmd, shell=True)

        shutil.rmtree('./tmp/')

    def analysis_species_subspecies(self, species_name, save_dir='.', save_align_png=True, make_phylogenetic_tree=False, bootstrap_times=100):
        self._make_tmp_dir()
        
        species_list = [species for species in self.taxonomy2otu['species'].keys() if species_name in species]
        
        for species in species_list: 
            otu_list = self.taxonomy2otu['species'][species]
            for otu in otu_list:
                seq = ''
                for uniq in self.otu2unique[otu]:
                    subseq = self.uniq_seq[uniq]
                    subseq = self._split_lines(subseq)
                    seq = seq + f'>{otu}\n{subseq}\n'
                with open(f'./tmp/{species}({otu}).fasta', 'w') as f:
                    f.write(seq)
        
        file_list = os.listdir('./tmp')
        file_list = [f'./tmp/{file}' for file in file_list]
        file_string = ' '.join(file_list)
        self._make_umap(file_string=file_string, save_folder_name=f'{save_dir}/{species_name}')

        seq = ''
        for species in species_list:
            otu_list = self.taxonomy2otu['species'][species]
            for otu in otu_list:
                subseq = self.otu_seq[otu]
                subseq = self._split_lines(subseq)
                seq = seq + f'>{otu}\n{subseq}\n'
        with open('./tmp/seq.fa', 'w') as f:
            f.write(seq)
        
        self._make_align_file(seq_file='./tmp/seq.fa', aln_file='./tmp/seq.aln')
        
        if os.path.exists('./tmp/seq.aln'):
            self._show_alignment(aln_seq_path='./tmp/seq.aln', save_path=f'{save_dir}/{species_name}/alignment.png', save_png=save_align_png)
        else:
            print(f'{species_name} contains 1 sequence, nothing to align')

        if make_phylogenetic_tree==True:
            cmd = f'iqtree2 -m GTR+F+I+G4 -s ./tmp/seq.aln -b {bootstrap_times} --prefix {species_name}'
            subprocess.run(cmd, shell=True)
            cmd = f'move {species_name}.* {save_dir}/{species_name}'
            subprocess.run(cmd, shell=True)

        shutil.rmtree('./tmp/')

    def usum_sample(self, save_dir='.'):
        self._make_tmp_dir()

        for species in self.taxonomy2otu['species'].keys():
            seq = ''
            otu_list = self.taxonomy2otu['species'][species]
            for otu in otu_list:
                subseq = self.otu_seq[otu]
                subseq = self._split_lines(subseq)
                seq = seq + f'>{otu}\n{subseq}\n'
            with open(f'./tmp/{species}.fasta', 'w') as file:
                file.write(seq)

        file_list = ['./tmp/' + species + '.fasta' for species in self.taxonomy2otu['species']]
        file_string = ' '.join(file_list)
        self._make_umap(file_string=file_string, save_folder_name=f'{save_dir}/sample')

        shutil.rmtree('./tmp/')

    def _get_level_size(self, level_dict):
        level_size = {}

        for key, otu_list in level_dict.items():
            size = 0
            for otu in otu_list:
                size = size + int(self.otu_size[otu])
            level_size[key] = size

        total_size = sum(level_size.values())
        for key in level_size.keys():
            level_size[key] = level_size[key]/total_size * 100

        return level_size

    def barplot_sample(self, level, save_dir=None):
        level_dict = self.taxonomy2otu[level]
        level_size = self._get_level_size(level_dict)

        plotdata = pd.DataFrame(level_size, index=[f'sample'])
        fig = px.bar(plotdata, barmode='stack', labels={'value': 'Percentage (%)'},color_discrete_sequence=px.colors.qualitative.Pastel)
        fig.update_layout(xaxis_title="Sample No.", yaxis_title="Percentage (%)", legend=dict(x=1.05, y=1, traceorder='normal', orientation='h'))
        fig.show()
        if type(save_dir)=='str':
            fig.write_html(f'{save_dir}/sample_{level}_bar_chart.html')

class Otu(OtuAnalysis):
    def _read_otu_report(self, otu_report_path):
        with open(otu_report_path, 'r') as file:
            for line in file.readlines():
                line_list = line.split(';')
                haploid, uniq_size = line_list[0], line_list[1].replace('size=','')
                if 'otu' in line_list[2]:
                    uniq_top = re.search(r'otu\d+', line[2]).group(0).replace('o','O')
                elif 'perfect\s' in line[2] or 'match' in line[2] or 'noisy' in line[2] or 'perfect_chimera' in line[2] :
                    uniq_top = re.search(r'top=[a-zA-Z]+\d+', line).group(0).replace('top=','')
                
                if uniq_top not in self.otu2unique:
                    self.otu2unique[uniq_top] = []
                self.otu2unique[uniq_top].append(line[0].replace('Uniq',''))
                self.uniq_size[haploid] = uniq_size

class Zotu(OtuAnalysis):
    def _read_otu_report(self, otu_report_path):
        with open(otu_report_path, 'r') as file:
            zotu_num = 1
            chi_num = 1
            for line in file.readlines():
                line_list = line.split(';')
                haploid= line_list[0] 
                if 'denoise' in line_list[2]:
                    if 'amp' in line_list[2]:
                        uniq_top = haploid
                    elif 'shifted' in line_list[2] or 'bad' in line_list[2]:
                        uniq_top = re.search(r'top=Uniq\d+', line_list[3]).group(0).replace('top=','')

                    if uniq_top not in self.otu2unique:
                        self.otu2unique[uniq_top] = []
                    self.otu2unique[uniq_top].append(haploid)
                    self.uniq_size[haploid] = line_list[1].replace('size=','')

                elif 'chfilter' in line_list[2]:
                    if 'zotu' in line_list[2]:
                        k_new = f'Zotu{zotu_num}'
                        self.otu2unique[k_new] = self.otu2unique.pop(haploid)
                        zotu_num += 1
                    elif 'chimera' in line_list[2]:
                        k_new = f'Chimera{chi_num}'
                        self.otu2unique[k_new] = self.otu2unique.pop(haploid)
                        chi_num += 1

class SumAllSample(OtuAnalysis):

    def __init__(self, read_dir):
        self.sample_dict = {}
        self.read_dir = read_dir

        filename_list = os.listdir(f'{read_dir}/4_derep')
        self.samplenum_list = [filename.replace('_derep.fasta','') for filename in filename_list if filename.endswith('_derep.fasta')]
        for sample_num in self.samplenum_list:
            self.sample_dict[sample_num] = Zotu(read_dir=read_dir, sample_num=sample_num)

    def _read_otu_report(report_path):
        return super()._read_otu_report()

    def barplot_all(self, level, save=True):
        sample_size = {}
        level2otu = {}
        for sample_num in self.samplenum_list:
            sample_size[sample_num] = {}
            level2otu[sample_num] = self.sample_dict[sample_num].taxonomy2otu[level]
            for key, otu_list in level2otu[sample_num].items():
                size = 0
                for otu in otu_list:
                    size = size + self.sample_dict[sample_num].otu_size[otu]
                sample_size[sample_num][key] = size
            total_size = sum(sample_size[sample_num].values())
            for key in sample_size[sample_num].keys():
                sample_size[sample_num][key] = sample_size[sample_num][key]/total_size * 100

        samplesize_list = [sample_size[sample_num] for sample_num in self.samplenum_list]
        all_keys = list(set().union(*samplesize_list))
        all_keys.sort()

        for sample_num in self.samplenum_list:
            sample_size[sample_num] = [sample_size[sample_num].get(key, 0) for key in all_keys]

        plotdata = pd.DataFrame(sample_size, index=all_keys)
        fig = px.bar(plotdata.transpose(), barmode='stack', labels={'value': 'Percentage (%)'},color_discrete_sequence=px.colors.qualitative.Pastel)
        fig.update_xaxes(tickmode='linear')
        fig.update_layout(xaxis_title="Sample No.", yaxis_title="Percentage (%)", legend=dict(x=1.05, y=1, traceorder='normal', orientation='h'))
        fig.show()
        if save==True:
            fig.write_html(f'{level}_bar_chart.html')

    def analysis_species(self, species_name, save_dir = '.', save_align_png=True):
        self._make_tmp_dir()

        align_seq = ''
        file_string = ''
        for sample_num in self.samplenum_list:
            species2otu_list = self.sample_dict[sample_num].taxonomy2otu['species'].get(species_name)
            if species2otu_list is not None:
                for otu in species2otu_list:

                    umap_seq = ''
                    for uniq in self.sample_dict[sample_num].otu2unique[otu]:
                        seq = self.sample_dict[sample_num].uniq_seq[uniq]
                        seq = self._split_lines(seq)
                        umap_seq = umap_seq + f'>{uniq}\n{seq}\n'
                    with open(f'./tmp/{sample_num}_{otu}.fasta', 'w') as file:
                        file.write(umap_seq)
                    file_string = file_string + f'./tmp/{sample_num}_{otu}.fasta '

                    title = f'>{sample_num}_{otu}'
                    seq = self.sample_dict[sample_num].otu_seq[otu]
                    align_seq = align_seq + f'{title}\n{seq}\n'
        
        with open('./tmp/seq.fa', 'w') as file:
            file.write(align_seq)

        self._make_umap(file_string=file_string, save_folder_name=f'{save_dir}/{species_name}')

        self._make_align_file(seq_file='./tmp/seq.fa', aln_file=f'{save_dir}/{species_name}/{species_name}')
        
        if os.path.exists(f'{save_dir}/{species_name}/{species_name}.aln'):
            self._show_alignment(aln_seq_path='./tmp/seq.aln', save_path=f'{save_dir}/{species_name}/alignment.png', save_png=save_align_png)
        else:
            print(f'{species_name} contains 1 sequence, nothing to align')

        # shutil.rmtree('./tmp/')

    def analysis_species_subspecies(self, species_name, save_dir='.', neighbors = 15, umap_min_dist = 0.1, alignment=True, dry=False):
        self._make_tmp_dir()
        
        align_seq = ''
        file_string = ''
        species_list = [species for sample_num in self.samplenum_list for species in self.sample_dict[sample_num].taxonomy2otu['species'].keys() if species_name in species]
        species_set = set(species_list)
        if dry == False:
            for species in species_set:
                umap_seq = ''
                for sample_num in self.samplenum_list:
                    otu_list = self.sample_dict[sample_num].taxonomy2otu['species'].get(species, [])
                    if otu_list != []:
                        for otu in otu_list:
                            title = f'>{species}_{sample_num}_{otu}'
                            seq  = self.sample_dict[sample_num].otu_seq[otu]
                            umap_seq = umap_seq + f'{title}\n{seq}\n'
                            align_seq = align_seq + f'{title}\n{seq}\n'

                with open(f'./tmp/{species}.fasta', 'w') as file:
                    file.write(umap_seq)
                file_string = file_string + f'./tmp/{species}.fasta '

            with open(f'./tmp/seq.fa', 'w') as file:
                file.write(align_seq)

            self._make_umap(file_string=file_string, save_folder_name=f'{save_dir}/{species_name}', neighbors = neighbors, umap_min_dist = umap_min_dist)

            if alignment:
                self._make_align_file(seq_file='./tmp/seq.fa', aln_file=f'{save_dir}/{species_name}/{species_name}')

                if os.path.exists(f'{save_dir}/{species_name}/{species_name}.aln'):
                    self._show_alignment(aln_seq_path=f'{save_dir}/{species_name}/{species_name}.aln', save_path=f'{save_dir}/{species_name}/alignment.png')
                else:
                    print(f'{species_name} contains 1 sequence, nothing to align')

            shutil.rmtree('./tmp/')

    def umap_genus(self, genus_name, save_dir='.', neighbors = 15, umap_min_dist = 0.1, dry=False):
        self._make_tmp_dir()
        
        align_seq = ''
        file_string = ''
        species_list = [species for sample_num in self.samplenum_list for species in self.sample_dict[sample_num].taxonomy2otu['species'].keys() if genus_name in species]
        species_set = set(species_list)
        print(species_set)
        if dry == False:
            for species in species_set:
                umap_seq = ''
                for sample_num in self.samplenum_list:
                    otu_list = self.sample_dict[sample_num].taxonomy2otu['species'].get(species, [])
                    if otu_list != []:
                        for otu in otu_list:
                            title = f'>{species}_{sample_num}_{otu}'
                            seq  = self.sample_dict[sample_num].otu_seq[otu]
                            umap_seq = umap_seq + f'{title}\n{seq}\n'
                            align_seq = align_seq + f'{title}\n{seq}\n'

                with open(f'./tmp/{species}.fasta', 'w') as file:
                    file.write(umap_seq)
                file_string = file_string + f'./tmp/{species}.fasta '

            with open(f'./tmp/seq.fa', 'w') as file:
                file.write(align_seq)

            self._make_umap(file_string=file_string, save_folder_name=f'{save_dir}/{genus_name}', neighbors = neighbors, umap_min_dist = umap_min_dist)

        shutil.rmtree('./tmp/')

    def umap_family(self, family_target, save_dir='.', neighbors = 15, umap_min_dist = 0.1, dry=False):
        self._make_tmp_dir()

        print(family_target)
        species_seq = {} # key is species_name, value is its sequence from all samples
        file_string = ''
        for sample_num in self.samplenum_list:
            for family, genus_list in self.sample_dict[sample_num].taxonomylevel['family'].items():
                if family_target == family:
                    #to genus level
                    for genus_target in genus_list:
                        for genus, species_list in self.sample_dict[sample_num].taxonomylevel['genus'].items():
                            if genus_target == genus:
                                #to species level
                                for species_target in species_list:
                                    if species_target not in species_seq:
                                        species_seq[species_target] = ""
                                    # to zotu level
                                    zotu_list = self.sample_dict[sample_num].taxonomy2otu['species'][species_target]
                                    for zotu in zotu_list:
                                        # to haplotype level
                                        haplotype_list = self.sample_dict[sample_num].otu2unique[zotu]
                                        for haplotype in haplotype_list:
                                            title = f'>{species_target}_{sample_num}_{zotu}_{haplotype}'
                                            seq  = self.sample_dict[sample_num].uniq_seq[haplotype]
                                            species_seq[species_target] += f'{title}\n{seq}\n'
        if dry == False:
            for species, seq in species_seq.items():
                with open(f'./tmp/{species}.fasta', 'w') as file:
                    file.write(seq)
                file_string = file_string + f'./tmp/{species}.fasta '

            self._make_umap(file_string=file_string, save_folder_name=f'{save_dir}/{family_target}', neighbors = neighbors, umap_min_dist = umap_min_dist)

        shutil.rmtree('./tmp/')

    def species_multiple_otu(self):
        species2otu = {}
        for sample_num in self.samplenum_list:
            species2otu[sample_num] = self.sample_dict[sample_num].taxonomy2otu['species']
        species2otu_list = [species2otu[sample_num] for sample_num in species2otu]
        all_keys = list(set().union(*species2otu_list))
        allsample_species = {}
        for key in all_keys:
            allsample_species[key] = []
        for sample_num in self.samplenum_list:
            for species, otu_list in species2otu[sample_num].items():
                allsample_species[species].extend(otu_list)
        
        species_num = 0
        species_list = []
        for species, otu_list in allsample_species.items():
            if len(otu_list)>1 and len(species.split('_'))==2:
                species_num += 1
                species_list.append(species)
        print(species_num)
        print(species_list)


if __name__ == '__main__':
    # sample1 = Zotu(read_dir='./data/keelung/2303', sample_num='2303-H02')
    # print(sample1.taxonomy2otu)
    # sample1.within_otu_align('Zotu5', save=False)
    # sample1.analysis_species_subspecies(species_name='Mugil_cephalus')
    # sample16.usum_sample()
    # sample1.barplot_sample(level='family', save_dir=None)
    read_path = './data/all_site'
    # read_path = './taoyuan'
    # read_path = './keelung/3month'
    a = SumAllSample(read_dir=read_path)
    # species_list = [list(a.sample_dict[sample_num].taxonomy2otu['species'].keys()) for sample_num in a.samplenum_list]
    # all_species = list(set().union(*species_list))
    # for species in all_species:
    #     if len(species.split('_')) == 2:
    #         a.analysis_species_subspecies(species_name=species, dry=False, alignment=False, save_dir='./umap_result')
    family_list = [list(a.sample_dict[sample_num].taxonomy2otu['family'].keys()) for sample_num in a.samplenum_list]
    all_family = list(set().union(*family_list))
    for family in all_family:
        a.umap_family(family_target=family, save_dir='./result/umap_family_haplotype', dry=False)
    # for species in all_species:
    #     if len(species.split('_')) == 2:
    #         a.analysis_species_subspecies(species_name=species, dry=False, alignment=False, save_dir='./umap_result')

    # a.barplot_all(level='species', save=True)
    # a.analysis_species('Mugil_cephalus')
    # a.analysis_species_subspecies(species_name='Planiliza_macrolepis', neighbors=10, umap_min_dist=0, dry=False, alignment=False)
    # a.umap_genus(genus_name='Mugil', neighbors = 20, umap_min_dist = 0, dry=False)
    # a.species_multiple_otu()

