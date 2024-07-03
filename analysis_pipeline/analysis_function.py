import os
import subprocess
import plotly.express as px
import shutil
from typing import Dict, List
import pandas as pd
import plotly.graph_objects as go
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from Bio import SeqIO

# TODO(SW): This module is a collection of functions, but probably should be split into several smaller modules. Each module should have a clear purpose.



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
    fig.update_layout(xaxis_title="Sample ID", yaxis_title="Percentage (%)", legend={"x":1.05, "y":1, "traceorder":'normal', "orientation":'h'})
    return fig

def create_dir(dir_path):
    try:
        os.makedirs(dir_path, exist_ok=True)
        print("> Creating directory: {}".format(dir_path))
    except FileExistsError:
        print("> Directory already exists: {}".format(dir_path))
    except Exception as e:
        print("> Error creating directory: {}".format(dir_path))
        print(e)

def remove_dir(dir_path):
    # TODO(SW): try-except block?
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

def dereplicate_fasta(seq_file, uniq_file, relabel, threads=12, sizeout=False):
    cmd = f'usearch -fastx_uniques {seq_file} -threads {threads} \
            -relabel {relabel}- -fastaout {uniq_file}'
    if sizeout:
        cmd += ' -sizeout'
    print("> Running USEARCH command: ", cmd)
    try:
        subprocess.run(cmd, shell=True)
    except Exception as e:
        print(e)


def write_fasta(units2fasta_dict, seq_file, dereplicate=False, sizeout=False):
    # TODO(SW): Move this to a different module.
    if dereplicate:
        fasta_list = []
        for unit_name, unit_seq in units2fasta_dict.items():
            with open(f'./tmp/{unit_name}.fa', 'w') as file:
                file.write(unit_seq)
            dereplicate_fasta(seq_file=f'./tmp/{unit_name}.fa', uniq_file=f'./tmp/{unit_name}_uniq.fa', relabel=unit_name, sizeout=sizeout)
            with open(f'./tmp/{unit_name}_uniq.fa', 'r') as file:
                fasta_list.append(file.read())
        with open(seq_file, 'w') as file:
            fasta_str = "".join(fasta_list)
            file.write(fasta_str)
    else:
        fasta_list = list(units2fasta_dict.values())
        with open(seq_file, 'w') as file:
            fasta_str = "".join(fasta_list)
            file.write(fasta_str)
    num_seq = fasta_str.count(">")
    print(f"\n> Written {num_seq} sequences to:  {seq_file}")
    return num_seq

def align_fasta(seq_file, aln_file):
    cmd = f'clustalo -i {seq_file} -o {aln_file} --force'
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
            '-redo': ' -redo',
            '--redo-tree': ' --redo-tree',
            '--undo': ' --undo',
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
    # if not os.path.isdir(save_subdir):
    os.makedirs(save_subdir, exist_ok=True)

    model = model or 'TEST'

    bootstrap_string = f'-b {bootstrap} ' if bootstrap is not None else ''
    threads = threads or 'AUTO'

    # TODO(SW): Technically shell=True is a security risk. It's convient, so it's used here. But you should be aware of this.
    # TODO(SW): Try to bulid it wint shell=False and shlex.split().
    cmd = f'iqtree2 -m {model} -s {seq_path} {bootstrap_string}--prefix {os.path.join(save_subdir, prefix)} -nt {threads}{checkpoint}'
    print("> Running IQTREE2 command: ", cmd)
    subprocess.run(cmd, shell=True)

def remove_row_by_unit_occurance(index, n):
    counts = index["unit"].value_counts()
    units_to_remove = counts[counts < n].index
    index = index[~index["unit"].isin(units_to_remove)]
    return index

def get_color_hex(n, cmap="rainbow"):
    color_key = plt.get_cmap(cmap)(np.linspace(0, 1, n))
    color_hex = [matplotlib.colors.to_hex(color_key[i]) for i in range(n)]
    print(f"> Color key: {" ".join(color_hex)}")
    return color_hex

def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def get_uniq_seq_freq(seq_file, uniq_seq_file, seq_labels):
    #TODO(SW): Needs to close() these files.
    uniq_labels = np.unique(seq_labels)
    unit_uniq_seq = SeqIO.parse(open(uniq_seq_file), 'fasta')
    label_freq_each_uniq_seq = {}
    for uniq_seq in unit_uniq_seq:
        label_freq_each_uniq_seq[uniq_seq.name] = {label:0 for label in uniq_labels}
        unit_seq = SeqIO.parse(open(seq_file), 'fasta')
        for i, seq in enumerate(unit_seq):
            if seq.seq == uniq_seq.seq:
                label_freq_each_uniq_seq[uniq_seq.name][seq_labels[i]] += 1

    freq_string = f"""
Begin Traits;
Dimensions NTraits={len(uniq_labels)};
format labels=yes missing=? separator=Comma;
TraitLabels {" ".join(map(str, uniq_labels))};
Matrix
"""
    for seq_id, freq_each_label in label_freq_each_uniq_seq.items():
        freq_string += f"{seq_id} {",".join(map(str, list(freq_each_label.values())))}\n"
    freq_string += """
;
end;
"""
    return freq_string