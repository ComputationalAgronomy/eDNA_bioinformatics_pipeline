from abc import ABC, abstractmethod
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re
import umap
from pymsaviz import MsaViz
# import matplotlib.colors as mcolors


class OtuAnalysis(ABC):

    def __init__(self):
        self.otu2uniq = {}
        self.uniq2otu = []
        self.otu2spc = {}
        self.uniq_seq = {}
        self.uniq_size = []
        self.uniq_type = []
        self.uniq_dqt = []
        self.otu_seq = {}
        self.otu_size = []

    def _read_seq(self, seq_path, prefix = '>'):
        with open(seq_path,'r') as file:
            seq={}
            for line in file.readlines():
                if prefix in line:
                    key_num = re.search(f'^{prefix}[a-zA-Z]+\d+', line).group(0)
                    key_num = re.sub(f'^{prefix}[a-zA-Z]+', '', key_num)
                    seq[key_num]=""
                else:
                    seq[key_num] = seq[key_num] + line.strip()
            return seq

    @abstractmethod
    def _read_otu_report(self, report_path):
        raise NotImplementedError

    def _get_otu_size(self):
        for uniq_group in self.otu2uniq.values():
            size = sum(int(self.uniq_size[int(uniq)-1]) for uniq in uniq_group)
            self.otu_size.append(size)

    # import derep.fasta, otu.fasta, otu_size.txt to get relationship between OTUs and unique sequences. 
    def import_data(self, uniq_path, otu_path, otu_report_path):
        self.uniq_seq = self._read_seq(uniq_path)
        self.otu_seq = self._read_seq(otu_path)
        self._read_otu_report(otu_report_path)

    # Split the sequences to conform to the fasta format.
    @staticmethod
    def _split_lines(seq):
        lines = [seq[i:i+59] for i in range(0, len(seq), 60)]
        return '\n'.join(lines)
    
    #plot MSA between unique sequences in one OTU.
    def within_otu_align(self, otu_name, save=False):
        uniq_list = self.otu2uniq[otu_name]
        seq = [self.uniq_seq[uniq] for uniq in uniq_list]
        text = ""
        for i, subseq in enumerate(seq):
            subseq = self._split_lines(subseq)
            text = f"{text}>Uniq{uniq_list[i]}\n{subseq}\n"
        with open('seq.fa', 'w') as f:
            f.write(text)
        
        cmd = '.\clustal-omega-1.2.2-win64\clustalo.exe -i seq.fa -o seq.aln'
        os.system(cmd)
        os.remove('seq.fa')
        mv = MsaViz('seq.aln', wrap_length=60, show_count=True, show_grid=False, show_consensus=True)
        mv.plotfig()
        if save==True:
            mv.savefig(f'{otu_name}.png')
        os.remove('seq.aln')

    # plot UMAP and alignment between OTUs in one species.
    def usum_otu(self, otu_list, spc_name):
        for otu_name in otu_list:
            seq = ''
            otu2uniq = self.otu2uniq[f'Zotu{otu_name}']
            for uniq in otu2uniq:
                subseq = self.uniq_seq[uniq]
                subseq = self._split_lines(subseq)
                seq = seq + f'>Uniq{uniq}\n{subseq}\n'
            with open(f'Zotu{otu_name}.fasta', 'w') as f:
                f.write(seq)
        cmd = 'usum ' 
        for otu in otu_list:
            cmd = cmd + f'Zotu{otu}.fasta '
        cmd = cmd +  '--labels ' 
        for otu in otu_list:
            cmd = cmd + f'Zotu{otu} '
        cmd = cmd + f'--maxdist 1.0 --termdist 1.0 --output {spc_name} -f'
        os.system(cmd)
        for otu in otu_list:
            os.remove(f'Zotu{otu}.fasta')
    
        seq=''
        for otu in otu_list:
            subseq = self.otu_seq[otu]
            subseq = self._split_lines(subseq)
            seq = seq + f'>Zotu{otu}\n\{subseq}\n'
        with open('seq.fa', 'w') as f:
            f.write(seq)
        cmd = '.\clustal-omega-1.2.2-win64\clustalo.exe -i seq.fa -o seq.aln'
        os.system(cmd)
        os.remove('seq.fa')
        mv = MsaViz('seq.aln', wrap_length=60, show_count=True, show_grid=False, show_consensus=True)
        mv.plotfig()
        mv.savefig(f'./{spc_name}/alignment.png')
        os.remove('seq.aln')

class Otu(OtuAnalysis):
    def _read_otu_report(self, otu_report_path):
        with open(otu_report_path, 'r') as file:
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
                if u_top not in self.otu2uniq:
                    self.otu2uniq[u_top] = []
                self.otu2uniq[u_top].append(line[0].replace('Uniq',''))
                self.uniq2otu.append(str(u_top))
                self.uniq_size.append(u_size)
                self.uniq_type.append(u_type)
            self._get_otu_size()

class Zotu(OtuAnalysis):
    def _read_otu_report(self, otu_report_path):
        with open(otu_report_path, 'r') as file:
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
                    if u_top not in self.otu2uniq:
                        self.otu2uniq[u_top] = []
                    self.otu2uniq[u_top].append(line[0].replace('Uniq',''))
                    self.uniq2otu.append(str(u_top))
                    self.uniq_size.append(u_size)
                    self.uniq_type.append(u_type)
                    self.uniq_dqt.append(dqt)
                else:
                    uniq_num = int(line[0].replace('Uniq',''))-1
                    if 'zotu' in line[2]:
                        k_new = f'Zotu{zotu_num}'
                        self.otu2uniq[k_new] = self.otu2uniq.pop(line[0])
                        self.uniq_type[uniq_num] = 'zotu'
                        zotu_num += 1
                    if 'chimera' in line[2]:
                        k_new = f'Chimera{chi_num}'
                        self.otu2uniq[k_new] = self.otu2uniq.pop(line[0])
                        self.uniq_type[uniq_num] = 'chimera'
                        self.uniq_dqt[uniq_num] = line[3].replace('dpt=','')
                        chi_num += 1
            self._get_otu_size()

if __name__ == '__main__':
    
    sample_num = 1
    uniq_path = f'./cleandata/4_derep/{sample_num}_derep.fasta'
    otu_path = f'./cleandata/5_haploid/otu/{sample_num}_otu.fasta'
    otu_report_path = f'./cleandata/5_haploid/otu/{sample_num}_otu_size.txt'
    zotu_path = f'./cleandata/5_haploid/zotu/{sample_num}_zotu.fasta'
    zotu_report_path = f'./cleandata/5_haploid/zotu/{sample_num}_zotu_size.txt'
    
    a = Zotu()
    a.import_data(uniq_path=uniq_path, otu_path=zotu_path, otu_report_path=zotu_report_path)
    
    otu_list = ['5', '17']
    species_name = 'Saurida_umeyoshii'
    a.usum_otu(otu_list=otu_list, spc_name=species_name)
    

