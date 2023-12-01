from abc import ABC, abstractmethod
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import subprocess
import re
from pymsaviz import MsaViz
import shutil
from matplotlib.cm import get_cmap
import plotly.express as px
import glob

class OtuAnalysis(ABC):

    def __init__(self, read_dir, sample_num):
        self.uniq_seq = {}
        self.uniq_size = {}
        self.otu_seq = {}
        self.otu_size = {}
        self.otu2unique = {}
        self.taxonomy2otu = {'kingdom':{}, 'phylum':{}, 'class':{}, 'order':{}, 'family':{}, 'genus':{}, 'species':{}}

        uniq_path = f'{read_dir}/4_derep/{sample_num}_derep.fasta'
        otu_path = f'{read_dir}/5_haploid/zotu/{sample_num}_zotu.fasta'
        otu_report_path = f'{read_dir}/5_haploid/zotu/{sample_num}_zotu_size.txt'
        blast_path = f'{read_dir}/6_blastn/mifish_db/{sample_num}_zotu.csv'
        
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
            self.taxonomy2otu['family']['no_family_level'] = []
            self.taxonomy2otu['order']['no_order_level'] = []

            for line in file.readlines():
                blast_list = line.split(',')
                haploid, species, genus, family, order, class_, phylum, kingdom = blast_list[0], blast_list[2], blast_list[3], blast_list[4], blast_list[5], blast_list[6], blast_list[7], blast_list[8]

                if species not in self.taxonomy2otu['species']:
                    self.taxonomy2otu['species'][species] = []
                self.taxonomy2otu['species'][species].append(haploid)

                if genus not in self.taxonomy2otu['genus']:
                    self.taxonomy2otu['genus'][genus] = []
                self.taxonomy2otu['genus'][genus].append(haploid)

                if family=='':
                    self.taxonomy2otu['family']['no_family_level'].append(haploid)
                else:
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

    # Split the sequences to conform to the fasta format.
    @staticmethod
    def _split_lines(seq):
        lines = [seq[i:i+59] for i in range(0, len(seq), 60)]
        return '\n'.join(lines)
    
    @staticmethod
    def make_tmp_dir():
        if not os.path.isdir('./tmp/'):
            os.makedirs('./tmp/')

    #plot MSA between unique sequences in one OTU.
    def within_otu_align(self, otu_name, save=False):
        self.make_tmp_dir()

        uniq_list = self.otu2unique[otu_name]
        text = ""
        for uniq in uniq_list:
            text = text + f'>{uniq}\n{self.uniq_seq[uniq]}\n'
        with open('./tmp/seq.fa', 'w') as f:
            f.write(text)
        
        cmd = 'clustalo -i ./tmp/seq.fa -o ./tmp/seq.aln'
        subprocess.Popen(cmd, shell=True)

        if os.path.exists('./tmp/seq.aln'):
            mv = MsaViz('./tmp/seq.aln', wrap_length=60, show_count=True, show_grid=False, show_consensus=True)
            mv.plotfig()
            if save==True:
                mv.savefig(f'{otu_name}.png')
        else:
            print(f'{otu_name} contains 1 sequence, nothing to align')
        
        shutil.rmtree('./tmp/')

    # plot UMAP and alignment between OTUs in one species.
    def analysis_species(self, species_name, make_phylogenetic_tree=False, bootstrap_times=100):
        self.make_tmp_dir()

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
        otu_string = ' '.join(otu_list)
        cmd = f'usum {file_string} --labels {otu_string} --maxdist 1.0 --termdist 1.0 --output {species_name} -f'
        subprocess.Popen(cmd, shell=True)

        seq = ''
        for otu in otu_list:
            subseq = self.otu_seq[otu]
            subseq = self._split_lines(subseq)
            seq = seq + f'>{otu}\n{subseq}\n'
        with open('./tmp/seq.fa', 'w') as f:
            f.write(seq)
        cmd = 'clustalo -i ./tmp/seq.fa -o ./tmp/seq.aln'
        subprocess.Popen(cmd, shell=True)

        mv = MsaViz('./tmp/seq.aln', wrap_length=60, show_count=True, show_grid=False, show_consensus=True)
        mv.plotfig()
        mv.savefig(f'./{species_name}/alignment.png')
        
        if make_phylogenetic_tree==True:
            cmd = f'iqtree2 -m GTR+F+I+G4 -s ./tmp/seq.aln -b {bootstrap_times} --prefix {species_name}'
            subprocess.Popen(cmd, shell=True)
            cmd = f'move {species_name}.* {species_name}'
            subprocess.Popen(cmd, shell=True)

        # shutil.rmtree('./tmp/')

    def analysis_species_subspecies(self, species_name, make_phylogenetic_tree=False, bootstrap_times=100):
        self.make_tmp_dir()
        
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
        cmd = f'usum {file_string} --maxdist 1.0 --termdist 1.0 --output {species_name} -f'
        subprocess.Popen(cmd, shell=True)

        seq = ''
        for species in species_list:
            otu_list = self.taxonomy2otu['species'][species]
            for otu in otu_list:
                subseq = self.otu_seq[otu]
                subseq = self._split_lines(subseq)
                seq = seq + f'>{otu}\n{subseq}\n'
        with open('./tmp/seq.fa', 'w') as f:
            f.write(seq)
        cmd = 'clustalo -i ./tmp/seq.fa -o ./tmp/seq.aln'
        subprocess.Popen(cmd, shell=True)
        mv = MsaViz('./tmp/seq.aln', wrap_length=60, show_count=True, show_grid=False, show_consensus=True)
        mv.plotfig()
        mv.savefig(f'./{species_name}/alignment.png')
        
        if make_phylogenetic_tree==True:
            cmd = f'iqtree2 -m GTR+F+I+G4 -s ./tmp/seq.aln -b {bootstrap_times} --prefix {species_name}'
            subprocess.Popen(cmd, shell=True)
            cmd = f'move {species_name}.* {species_name}'
            subprocess.Popen(cmd, shell=True)

        shutil.rmtree('./tmp/')

    def usum_sample(self):
        self.make_tmp_dir()

        for species in self.taxonomy2otu['species']:
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
        species_string = ' '.join(self.taxonomy2otu['species'].keys())
        cmd = f'usum {file_string} --labels {species_string} --maxdist 1.0 --termdist 1.0 --output sample -f'
        subprocess.Popen(cmd, shell=True)
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

    def barplot_sample(self, level, save=True):
        level_dict = self.taxonomy2otu[level]
        level_size = self._get_level_size(level_dict)

        plotdata = pd.DataFrame(level_size, index=[f'sample'])
        fig = px.bar(plotdata, barmode='stack', labels={'value': 'Percentage (%)'},color_discrete_sequence=px.colors.qualitative.Pastel)
        fig.update_layout(xaxis_title="Sample No.", yaxis_title="Percentage (%)", legend=dict(x=1.05, y=1, traceorder='normal', orientation='h'))
        fig.show()
        if save==True:
            fig.write_html(f'sample_{level}_bar_chart.html')

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

class SumAllSample:
    def __init__(self):
        self.sample_dict = {}

    def add_sample(self, read_dir, samplenum_list):
        for sample_num in samplenum_list:
            self.sample_dict[sample_num] = Zotu(read_dir='./cleandata', sample_num=sample_num)

    def barplot_all(self, level, samplenum_list, save=True):
        sample_size = {}
        level2otu = {}
        for sample_num in samplenum_list:
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

        samplesize_list = [sample_size[sample_num] for sample_num in samplenum_list]
        all_keys = list(set().union(*samplesize_list))
        all_keys.sort()

        for sample_num in samplenum_list:
            sample_size[sample_num] = [sample_size[sample_num].get(key, 0) for key in all_keys]

        plotdata = pd.DataFrame(sample_size, index=all_keys)
        fig = px.bar(plotdata.transpose(), barmode='stack', labels={'value': 'Percentage (%)'},color_discrete_sequence=px.colors.qualitative.Pastel)
        fig.update_xaxes(tickmode='linear')
        fig.update_layout(xaxis_title="Sample No.", yaxis_title="Percentage (%)", legend=dict(x=1.05, y=1, traceorder='normal', orientation='h'))
        fig.show()
        if save==True:
            fig.write_html(f'{level}_bar_chart.html')

if __name__ == '__main__':
    sample16 = Zotu(read_dir='./cleandata', sample_num=16)
    # sample1.within_otu_align('Zotu5', save=False)
    sample16.analysis_species_subspecies(species_name='Carassius_auratus')
    # sample1.usum_sample()
    # sample1.barplot_sample(level='family', save=False)
    # sample1.barplot_all('family', save=True)
    # a = SumAllSample()
    # a.add_sample(read_dir='./cleandata', samplenum_list=range(1,19))
    # a.barplot_all(level='species', samplenum_list=range(1, 19), save=False)
