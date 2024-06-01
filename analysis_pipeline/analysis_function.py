import os
import subprocess
import plotly.express as px
import shutil

def normalize_abundance(abundance_dict):
    total_size = sum(abundance_dict.values())
    norm_abundance = {key: value/total_size * 100 for key, value in abundance_dict.items()}
    return norm_abundance

def list_union(input_list): # e.g. [[1,2,3], [2,3,4], [3,4,5]] -> [1,2,3,4,5]
    uniq_list = list(set().union(*input_list))
    uniq_list.sort()
    return uniq_list

def create_barchart_fig(data):
    fig = px.bar(data, barmode='stack', labels={'value': 'Percentage (%)'},color_discrete_sequence=px.colors.qualitative.Pastel)
    fig.update_xaxes(tickmode='linear')
    fig.update_layout(xaxis_title="Sample ID", yaxis_title="Percentage (%)", legend=dict(x=1.05, y=1, traceorder='normal', orientation='h'))
    return fig

def create_dir(dir_path):
    if not os.path.isdir(dir_path):
        print(f"> Creating directory: {dir_path}")
        os.makedirs(dir_path)

def remove_dir(dir_path):
    if os.path.isdir(dir_path):
        print(f"> Removing directory: {dir_path}")
        shutil.rmtree(dir_path)

def write_umap_fasta(units2seq_dict):
    fasta_path_string = ""
    for unit_name, seq in units2seq_dict.items():
        with open(f'./tmp/{unit_name}.fa', 'w') as file:
            file.write(seq)
        fasta_path_string += f'./tmp/{unit_name}.fa '
    print(f"\n> Written fasta files to: {fasta_path_string}")
    return fasta_path_string # e.g. "./tmp/SpA.fa ./tmp/SpB.fa ./tmp/SpC.fa "

def umap_command(fasta_path_string, save_dir, prefix, neighbors, min_dist):
    save_path = os.path.join(save_dir, prefix)
    cmd = f'usum {fasta_path_string} --neighbors {neighbors} --umap-min-dist {min_dist} --maxdist 1.0 --termdist 1.0 --output {save_path} -f --seed 1'
    subprocess.run(cmd, shell=True)

def dereplicate_fasta(seq_file, uniq_file, relabel, threads=12):
    cmd = f'usearch -fastx_uniques {seq_file} -threads {threads} \
            -relabel {relabel}_ -fastaout {uniq_file}'
    print("> Running USEARCH command: ", cmd)
    subprocess.run(cmd, shell=True)

def write_mltree_fasta(units2fasta_dict, dereplicate=False):
    fasta_string = ""
    if dereplicate == True:
        for unit_name, unit_seq in units2fasta_dict.items():
            with open(f'./tmp/{unit_name}.fa', 'w') as file:
                file.write(unit_seq)
            dereplicate_fasta(seq_file=f'./tmp/{unit_name}.fa', uniq_file=f'./tmp/{unit_name}_uniq.fa', relabel=unit_name)
            with open(f'./tmp/{unit_name}_uniq.fa', 'r') as file:
                fasta_string += file.read()
        with open(f'./tmp/mltree.fa', 'w') as file:
            file.write(fasta_string)
    else:    
        for unit_seq in units2fasta_dict.values():
            fasta_string += unit_seq
        with open(f'./tmp/mltree.fa', 'w') as file:
            file.write(fasta_string)
    print(f"\n> Written MLTree fasta file to:  ./tmp/mltree.fa")

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
        while True:
            user_input = input("(-redo/--redo-tree/--undo/stop): ")
            if user_input in ['-redo', '--redo-tree', '--undo', 'stop', 'Stop', 'STOP']:
                return user_input
            else:
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
