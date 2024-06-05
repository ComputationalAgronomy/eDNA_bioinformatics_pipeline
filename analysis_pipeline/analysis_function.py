import os
import subprocess
import plotly.express as px
import shutil
from typing import Dict, List
import pandas as pd
import plotly.graph_objects as go

def normalize_abundance(abundance_dict: Dict[str, float]) -> Dict[str, float]:
    total_size = sum(abundance_dict.values())
    norm_abundance = {key: value/total_size * 100 for key, value in abundance_dict.items()}
    return norm_abundance

def list_union(lists_to_union: List[List[int]]) -> List[int]:
    uniq_list = list(set().union(*lists_to_union))
    uniq_list.sort()
    return uniq_list

def create_barchart_fig(data: pd.DataFrame) -> go.Figure:
    fig = px.bar(data, barmode='stack', labels={'value': 'Percentage (%)'},color_discrete_sequence=px.colors.qualitative.Pastel)
    fig.update_xaxes(tickmode='linear')
    fig.update_layout(xaxis_title="Sample ID", yaxis_title="Percentage (%)", legend=dict(x=1.05, y=1, traceorder='normal', orientation='h'))
    return fig

def create_dir(dir_path):
    if not os.path.isdir(dir_path):
        print("> Creating directory: {}"
              .format(dir_path))
        os.makedirs(dir_path)

def remove_dir(dir_path):
    if os.path.isdir(dir_path):
        print("> Removing directory: {}"
              .format(dir_path))
        shutil.rmtree(dir_path)

# def write_umap_fasta(units2seq_dict):
#     fasta_path_string = ""
#     for unit_name, seq in units2seq_dict.items():
#         with open(f'./tmp/{unit_name}.fa', 'w') as file:
#             file.write(seq)
#         fasta_path_string += f'./tmp/{unit_name}.fa '
#     print(f"\n> Written fasta files to: {fasta_path_string}")
#     return fasta_path_string # e.g. "./tmp/SpA.fa ./tmp/SpB.fa ./tmp/SpC.fa "

# def umap_command(fasta_path_string, save_dir, prefix, neighbors, min_dist):
#     save_path = os.path.join(save_dir, prefix)
#     cmd = f'usum {fasta_path_string} --neighbors {neighbors} --umap-min-dist {min_dist} --maxdist 1.0 --termdist 1.0 --output {save_path} -f --seed 1'
#     subprocess.run(cmd, shell=True)

def dereplicate_fasta(seq_file, uniq_file, relabel, threads=12):
    cmd = f'usearch -fastx_uniques {seq_file} -threads {threads} \
            -relabel {relabel}- -fastaout {uniq_file}'
    print("> Running USEARCH command: ", cmd)
    subprocess.run(cmd, shell=True)

def write_fasta(units2fasta_dict, seq_file, dereplicate=False):
    if dereplicate:
        fasta_strings = []
        for unit_name, unit_seq in units2fasta_dict.items():
            with open(f'./tmp/{unit_name}.fa', 'w') as file:
                file.write(unit_seq)
            dereplicate_fasta(seq_file=f'./tmp/{unit_name}.fa', uniq_file=f'./tmp/{unit_name}_uniq.fa', relabel=unit_name)
            with open(f'./tmp/{unit_name}_uniq.fa', 'r') as file:
                fasta_strings.append(file.read())
        with open(seq_file, 'w') as file:
            file.write("".join(fasta_strings))
    else:    
        fasta_string = "".join(units2fasta_dict.values())
        with open(seq_file, 'w') as file:
            file.write(fasta_string)
    print(f"\n> Written fasta file to:  {seq_file}")

def align_fasta(seq_file, aln_file):
    cmd = f'clustalo -i {seq_file} -o {aln_file} --distmat-out={aln_file}.mat --guidetree-out={aln_file}.dnd --full --force'
    print("Running Clustal Omega command: ", cmd)
    subprocess.run(cmd, shell=True)
    print(f"\n> Aligned fasta file to: {aln_file}\n")

def check_mltree_overwrite(save_dir, prefix):
    ckp_path = os.path.join(save_dir, prefix, prefix + '.ckp.gz') # e.g. save/dir/SpA/SpA.ckp.gz
    if os.path.exists(ckp_path):
        print(f"""
> MLTree checkpoint fileCheckpoint ({ckp_path}) indicates that a previous run successfully finished  already exists."
Use `-redo` option if you really want to redo the analysis and overwrite all output files.
Use `--redo-tree` option if you want to restore ModelFinder and only redo tree search.
Use `--undo` option if you want to continue previous run when changing/adding options.
        """)
        user_input_map = {
            '-redo': '-redo',
            '--redo-tree': '--redo-tree',
            '--undo': '--undo',
            'stop': 'stop',
            'Stop': 'stop',
            'STOP': 'stop'
        }
        while True:
            user_input = input("(-redo/--redo-tree/--undo/stop): ")
            if user_input in user_input_map:
                return user_input_map[user_input]
            print("> Invalid input.")
    else:
        return ""

def iqtree2_command(seq_path, save_dir, prefix, model, bootstrap, threads, checkpoint):
    save_subdir = os.path.join(save_dir, prefix)
    if not os.path.isdir(save_subdir):
        os.makedirs(save_subdir)
    model = model if model != None else 'TEST'
    bootstrap_string = f'-b {bootstrap} ' if bootstrap != None else ''
    threads = threads if threads != None else 'AUTO'
    cmd = f'iqtree2 -m {model} -s {seq_path} {bootstrap_string}--prefix {os.path.join(save_subdir, prefix)} -nt {threads}{checkpoint}'
    print("> Running IQTREE2 command: ", cmd)
    subprocess.run(cmd, shell=True)

def remove_row_by_unit_occurance(index, n):
    counts = index["unit"].value_counts()
    units_to_remove = counts[counts < n].index
    index = index[~index["unit"].isin(units_to_remove)]
    return index