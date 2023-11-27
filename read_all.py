from abc import ABC, abstractmethod
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import re
from pymsaviz import MsaViz
import shutil
from matplotlib.cm import get_cmap
import plotly.express as px

class OtuAnalysis(ABC):

    def __init__(self):
        self.sample = {}
        # self.otu2unique = {}
        # self.taxonomy = {}
        # self.species2otu = {}
        # self.genus2otu = {}
        # self.family2otu = {}
        # self.order2otu = {}
        # self.uniq_seq = {}
        # self.uniq_size = []
        # self.uniq_type = []
        # self.uniq_dqt = []
        # self.otu_seq = {}
        # self.otu_size = []

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

    @staticmethod
    @abstractmethod
    def _read_otu_report(report_path):
        raise NotImplementedError

    def _get_otu_size(self, num):
        otu_size = []
        for uniq_group in self.sample[num]['otu2unique'].values():
            size = sum(int(self.sample[num]['uniq_size'][int(uniq)-1]) for uniq in uniq_group)
            otu_size.append(size)
        return otu_size

    def _read_blast(self, blast_path):
        with open(blast_path, 'r') as file:
            species2otu, genus2otu, family2otu, order2otu, class2otu, phylum2otu, kingdom2otu = {}, {}, {}, {}, {}, {}, {}
            family2otu['no_family_level'] = []
            for line in file.readlines():
                blast_list = line.split(',')
                haploid, species, genus, family, order, class_, phylum, kingdom = blast_list[0], blast_list[2], blast_list[3], blast_list[4], blast_list[5], blast_list[6], blast_list[7], blast_list[8]
                if species not in species2otu:
                    species2otu[species] = []
                species2otu[species].append(haploid)

                if genus not in genus2otu:
                    genus2otu[genus] = []
                genus2otu[genus].append(haploid)

                if family=='':
                    family2otu['no_family_level'].append(haploid)
                else:
                    if family not in family2otu:
                        family2otu[family] = []
                    family2otu[family].append(haploid)

                if order not in order2otu:
                    order2otu[order] = []
                order2otu[order].append(haploid)

                if class_ not in class2otu:
                    class2otu[class_] = []
                class2otu[class_].append(haploid)

                if phylum not in phylum2otu:
                    phylum2otu[phylum] = []
                phylum2otu[phylum].append(haploid)

                if kingdom not in kingdom2otu:
                    kingdom2otu[kingdom] = []
                kingdom2otu[kingdom].append(haploid)

        return species2otu, genus2otu, family2otu, order2otu, class2otu, phylum2otu, kingdom2otu

    # import derep.fasta, otu.fasta, otu_size.txt to get relationship between OTUs and unique sequences. 
    def import_data(self, read_dir):
        for num in range(1, 19):
            uniq_path = f'{read_dir}/4_derep/{num}_derep.fasta'
            otu_path = f'{read_dir}/5_haploid/zotu/{num}_zotu.fasta'
            otu_report_path = f'{read_dir}/5_haploid/zotu/{num}_zotu_size.txt'
            blast_path = f'{read_dir}/6_blastn/mifish_db/{num}_zotu.csv'

            uniq_seq = self._read_seq(uniq_path)
            otu_seq = self._read_seq(otu_path)
            otu2unique, uniq_size, uniq_type, uniq_dqt = self._read_otu_report(otu_report_path)
            self.sample[num]={'uniq_seq':uniq_seq, 'otu_seq':otu_seq, 'otu2unique':otu2unique, 'uniq_size':uniq_size, 'uniq_type':uniq_type, 'uniq_dqt':uniq_dqt}
            self.sample[num]['otu_size'] = self._get_otu_size(num)
            species2otu, genus2otu, family2otu, order2otu, class2otu, phylum2otu, kingdom2otu= self._read_blast(blast_path)
            level_dict = {'species2otu':species2otu, 'genus2otu':genus2otu, 'family2otu':family2otu, 'order2otu':order2otu, 'class2otu':class2otu, 'phylum2otu':phylum2otu, 'kingdom2otu':kingdom2otu}
            self.sample[num].update(level_dict)
    
    # Split the sequences to conform to the fasta format.
    @staticmethod
    def _split_lines(seq):
        lines = [seq[i:i+59] for i in range(0, len(seq), 60)]
        return '\n'.join(lines)
    
    #plot MSA between unique sequences in one OTU.
    def within_otu_align(self, sample_num, otu_name, save=False):
        uniq_list = self.sample[sample_num]['otu2unique'][otu_name]
        seq = [self.sample[sample_num]['uniq_seq'][f'Uniq{uniq}'] for uniq in uniq_list]
        text = ""
        for i, subseq in enumerate(seq):
            subseq = self._split_lines(subseq)
            text = f"{text}>Uniq{uniq_list[i]}\n{subseq}\n"
        with open('seq.fa', 'w') as f:
            f.write(text)
        
        cmd = 'clustalo -i seq.fa -o seq.aln'
        os.system(cmd)
        os.remove('seq.fa')
        if os.path.exists('seq.aln'):
            mv = MsaViz('seq.aln', wrap_length=60, show_count=True, show_grid=False, show_consensus=True)
            mv.plotfig()
            if save==True:
                mv.savefig(f'{otu_name}.png')
            os.remove('seq.aln')
        else:
            print(f'{otu_name} contains 1 sequence, nothing to align')

    # plot UMAP and alignment between OTUs in one species.
    def analysis_species(self, sample_num, species_name):
        if not os.path.isdir('./tmp/'):
            os.makedirs('./tmp/')
        otu_list = self.sample[sample_num]['species2otu'][species_name]
        for otu_name in otu_list:
            seq = ''
            otu2unique = self.sample[sample_num]['otu2unique'][otu_name]
            for uniq in otu2unique:
                subseq = self.sample[sample_num]['uniq_seq'][f'Uniq{uniq}']
                subseq = self._split_lines(subseq)
                seq = seq + f'>Uniq{uniq}\n{subseq}\n'
            with open(f'./tmp/{otu_name}.fasta', 'w') as f:
                f.write(seq)
        file_list = ['./tmp/' + otu + '.fasta' for otu in otu_list]
        file_string = ' '.join(file_list)
        otu_string = ' '.join(otu_list)
        cmd = f'usum {file_string} --labels {otu_string} --maxdist 1.0 --termdist 1.0 --output {species_name} -f'
        os.system(cmd)
    
        seq = ''
        for otu in otu_list:
            subseq = self.sample[sample_num]['otu_seq'][otu]
            subseq = self._split_lines(subseq)
            seq = seq + f'>{otu}\n{subseq}\n'
        with open('./tmp/seq.fa', 'w') as f:
            f.write(seq)
        cmd = 'clustalo -i ./tmp/seq.fa -o ./tmp/seq.aln'
        os.system(cmd)
        mv = MsaViz('./tmp/seq.aln', wrap_length=60, show_count=True, show_grid=False, show_consensus=True)
        mv.plotfig()
        mv.savefig(f'./{species_name}/alignment.png')
        
        cmd = f'iqtree2 -m GTR+F+I+G4 -s ./tmp/seq.aln -b 1000 --prefix {species_name}'
        os.system(cmd)
        cmd = f'move {species_name}.* {species_name}'
        os.system(cmd)

        shutil.rmtree('./tmp/')

    def usum_sample(self, sample_num):
        if not os.path.isdir('./tmp/'):
            os.makedirs('./tmp/')
        for spc_name in self.sample[sample_num]['species2otu']:
            seq = ''
            otu_list = self.species2otu[spc_name]
            for otu in otu_list:
                subseq = self.sample[sample_num]['otu_seq'][otu]
                subseq = self._split_lines(subseq)
                seq = seq + f'>{otu}\n{subseq}\n'
            with open(f'./tmp/{spc_name}.fasta', 'w') as f:
                f.write(seq)

        file_list = ['./tmp/' + spc + '.fasta' for spc in self.sample[sample_num]['species2otu']]
        file_string = ' '.join(file_list)
        spc_string = ' '.join(self.sample[sample_num]['species2otu'].keys())
        cmd = f'usum {file_string} --labels {spc_string} --maxdist 1.0 --termdist 1.0 --output sample{sample_num} -f'
        os.system(cmd)
        shutil.rmtree('./tmp/')

    def _get_level_size(self, sample_num, level_dict):
        level_size = {}
        for key, otu_list in level_dict.items():
            size = 0
            for otu in otu_list:
                x = int(re.findall(r'\d+', otu)[0])-1
                size = size + self.sample[sample_num]['otu_size'][x]
            level_size[key] = size
        total_size = sum(level_size.values())
        for key in level_size.keys():
            level_size[key] = level_size[key]/total_size * 100
        return level_size

    def barplot_sample(self, sample_num, level, save=True):
        level_dict = self.sample[sample_num][f'{level}2otu']
        level_size = self._get_level_size(sample_num, level_dict)

        plotdata = pd.DataFrame(level_size, index=[f'sample{sample_num}'])
        fig = px.bar(plotdata, barmode='stack', labels={'value': 'Percentage (%)'},color_discrete_sequence=px.colors.qualitative.Pastel)
        fig.update_layout(xaxis_title="Sample No.", yaxis_title="Percentage (%)", legend=dict(x=1.05, y=1, traceorder='normal', orientation='h'))
        fig.show()
        if save==True:
            fig.write_html(f'{level}_sample{sample_num}_bar_chart.html')
        
    def barplot_all(self, level, save=True):
        sample_dict = {}
        for num in range(1, 19):
            level_dict = self.sample[num][f'{level}2otu']
            level_size = self._get_level_size(sample_num=num, level_dict=level_dict)
            sample_dict[f'sample{num}'] = level_size
        
        sample_list = [sample_dict[f'sample{num}'] for num in range(1,19)]
        all_key = list(set().union(*sample_list))
        all_key.sort()

        for num in range(1, 19):
            sample_dict[f'sample{num}'] = [sample_dict[f'sample{num}'].get(key, 0) for key in all_key]

        plotdata = pd.DataFrame(sample_dict, index=all_key)
        fig = px.bar(plotdata.transpose(), barmode='stack', labels={'value': 'Percentage (%)'},color_discrete_sequence=px.colors.qualitative.Pastel)
        fig.update_layout(xaxis_title="Sample No.", yaxis_title="Percentage (%)", legend=dict(x=1.05, y=1, traceorder='normal', orientation='h'))
        fig.show()
        if save==True:
            fig.write_html(f'{level}_bar_chart.html')

class Otu(OtuAnalysis):
    def _read_otu_report(self, otu_report_path):
        with open(otu_report_path, 'r') as file:
            otu2unique, uniq_size, uniq_type, uniq_dqt = {}, [], [], []
            for line in file.readlines():
                line = line.split(';')
                u_size = line[1].replace('size=','')
                if 'otu' in line[2]:
                    u_type = 'otu'
                    dqt = 0
                    u_top = re.search(r'otu\d+', line[2]).group(0).replace('o','O')
                if 'perfect\s' in line[2]:
                    u_type = 'perfect'
                    dqt = 0
                    u_top = re.search(r'top=[a-zA-Z]+\d+', line[2]).group(0).replace('top=','')
                if 'match' in line[2] or 'noisy' in line[2] or 'perfect_chimera' in line[2]:
                    u_type = 'match' if 'match' in line[2] else u_type
                    u_type = 'noisy_chimera' if 'noisy' in line[2] else u_type
                    u_type = 'perfect_chimera' if 'perfect_chimera' in line[2] else u_type
                    dqt = re.search(r'dqt=\d+', line[2]).group(0).replace('dqt=','')
                    top = re.search(r'top=[a-zA-Z]+\d+', line[3])
                    u_top = top.group(0).replace('top=','') if top is not None else -99
                if u_top not in otu2unique:
                    otu2unique[u_top] = []
                otu2unique[u_top].append(line[0].replace('Uniq',''))
                uniq_size.append(u_size)
                uniq_type.append(u_type)
                uniq_dqt.append(dqt)
        return otu2unique, uniq_size, uniq_type, uniq_dqt

class Zotu(OtuAnalysis):
    @staticmethod
    def _read_otu_report(otu_report_path):
        with open(otu_report_path, 'r') as file:
            otu2unique, uniq_size, uniq_type, uniq_dqt = {}, [], [], []
            zotu_num = 1
            chi_num = 1
            for line in file.readlines():
                line = line.split(';')
                u_size = line[1].replace('size=','')
                if 'amp' in line[2]:
                    u_type = ''
                    dqt = 0
                    u_top = line[0]
                if 'shifted' in line[2]:
                    u_type = 'shifted'
                    dqt = 0
                    u_top = re.search(r'top=Uniq\d+', line[3]).group(0).replace('top=','')
                if 'bad' in line[2]:
                    u_type = 'bad'
                    dqt = re.search(r'dqt=\d+', line[2]).group(0).replace('dqt=','')
                    u_top = re.search(r'top=Uniq\d+', line[3]).group(0).replace('top=','')
                if 'chfilter' not in line[2]:
                    if u_top not in otu2unique:
                        otu2unique[u_top] = []
                    otu2unique[u_top].append(line[0].replace('Uniq',''))
                    uniq_size.append(u_size)
                    uniq_type.append(u_type)
                    uniq_dqt.append(dqt)
                else:
                    uniq_num = int(line[0].replace('Uniq',''))-1
                    if 'zotu' in line[2]:
                        k_new = f'Zotu{zotu_num}'
                        otu2unique[k_new] = otu2unique.pop(line[0])
                        uniq_type[uniq_num] = 'zotu'
                        zotu_num += 1
                    if 'chimera' in line[2]:
                        k_new = f'Chimera{chi_num}'
                        otu2unique[k_new] = otu2unique.pop(line[0])
                        uniq_type[uniq_num] = 'chimera'
                        uniq_dqt[uniq_num] = line[3].replace('dpt=','')
                        chi_num += 1
        return otu2unique, uniq_size, uniq_type, uniq_dqt

if __name__ == '__main__':
    
    a = Zotu()
    a.import_data(read_dir='./cleandata')
    # a.within_otu_align(1, 'Zotu5', save=False)
    a.analysis_species(sample_num=1, species_name='Sardinella_fijiensis')
    # a.usum_sample()
    # a.barplot_sample(sample_num=2, level='family', save=False)
    # a.barplot_all('family', save=True)