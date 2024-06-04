from analysis_function import align_fasta, create_dir
from color_umap import plot_points
from run_hdbscan import *
import os
import numpy as np
import pandas as pd
from scipy import sparse
import subprocess
from Bio import SeqIO
import umap
import bokeh

def fasta2index(seq_file, fasta_file):
    index_list = []
    with open(fasta_file, 'w')as file:
        records = SeqIO.parse(seq_file, 'fasta')
        for i, record in enumerate(records):
            index = str(i)
            unit = record.description.split("-")[0]
            seq_id = record.description
            index_list.append([index, seq_id, unit])
            record.id = str(i)
            record.description = ''
            record.name = record.id
            SeqIO.write(record, file, 'fasta')
    
    index = pd.DataFrame(index_list, columns=["index", "seq_id", "unit"])
    return index

def calc_dist(seq_file, aln_file, dist_file, maxdist=1.0, termdist=1.0, threads=12):
    align_fasta(seq_file=seq_file, aln_file=aln_file)
    cmd = f"usearch -calc_distmx {aln_file} -tabbedout {dist_file} -maxdist {str(maxdist)} -termdist {str(termdist)}"
    print("> Running USEARCH command: ", cmd)
    if threads:
        cmd += f" -threads {str(threads)}"
    subprocess.run(cmd, shell=True)

def load_sparse_dist_matrix(dist_path):
    dist_matrix = pd.read_csv(dist_path, header=None, sep='\t')
    print(f'> Created sparse {max(dist_matrix[0])+1} x {max(dist_matrix[0])+1} distance matrix...')
    
    diagonal = dist_matrix[0] == dist_matrix[1]
    row = np.concatenate([dist_matrix[0], dist_matrix[1][~diagonal]])
    col = np.concatenate([dist_matrix[1], dist_matrix[0][~diagonal]])
    data = 1 - np.concatenate([dist_matrix[2], dist_matrix[2][~diagonal]])

    dist_matrix = sparse.csr_matrix((data, (row, col)), dtype=np.float32)
    return 1 - dist_matrix.toarray()

def fit_umap(dist_matrix, random_state=None, neighbors=15, min_dist=0.1, spread=1.0):
    print(f'> Creating UMAP embedding with {neighbors} neighbors...')
    reducer = umap.UMAP(n_neighbors=neighbors, random_state=random_state, min_dist=min_dist, spread=spread, metric='precomputed')
    embedding = reducer.fit_transform(dist_matrix)
    return reducer, embedding

def get_index_source_label(seq_id):
    s = []
    for id in seq_id:
        if "taoyuan" in id:
            s.append("taoyuan")
        elif "keelung" in id:
            s.append("keelung")
        else:
            s.append("reference")
    return s

def get_index_target_label(index, target2units):
    t = []
    for unit in index['unit']:
        for target, units in target2units.items():
            if unit in units:
                t.append(target)
    return t

def remove_row_by_unit_occurance(index, n):
    counts = index["unit"].value_counts()
    units_to_remove = counts[counts < n].index
    index = index[~index["unit"].isin(units_to_remove)]
    return index

def write_umap_file(seq_file, save_dir, target2units, random_state=42, neighbors=15, min_dist=0.1):
    create_dir(save_dir)
    index_fasta_file = os.path.join(save_dir, "input.fa")
    aln_file = os.path.join(save_dir, f"input.aln")
    dist_file = os.path.join(save_dir, "distance.txt")
    index_path = os.path.join(save_dir, "index.tsv")

    index = fasta2index(seq_file=seq_file, fasta_file=index_fasta_file)
    index['target'] = get_index_target_label(index, target2units)
    index['source'] = get_index_source_label(index['seq_id'])

    calc_dist(seq_file=index_fasta_file, aln_file=aln_file, dist_file=dist_file)
    dist_matrix = load_sparse_dist_matrix(dist_file)
    reducer, embedding = fit_umap(dist_matrix, random_state=random_state, neighbors=neighbors, min_dist=min_dist)
    index['umap1'] = embedding[:,0]
    index['umap2'] = embedding[:,1]

    index.to_csv(index_path, sep='\t', index=False)
    print(f'Saved index TSV to: {index_path}')

def plot_umap(index_file, save_dir, n_unit_threshold=1, theme='fire', width=800, height=800):
    png_path = os.path.join(save_dir, "umap.png") 
    html_path = os.path.join(save_dir, "umap.html")
    index = pd.read_csv(index_file, sep='\t')

    index = remove_row_by_unit_occurance(index, n=n_unit_threshold)

    points = index[["umap1", "umap2"]].to_numpy()

    print('\n> Drawing PNG...')
    ax = plot_points(points, labels=index['unit'], markers=index['source'], theme=theme, width=width, height=height)
    ax.figure.savefig(png_path, bbox_inches='tight')
    print(f'Saved PNG to: {png_path}')
    
    # print('\n> Drawing interactive plot...')
    # p = umap.plot.interactive(reducer, labels=index['label'], theme=theme, width=width, height=height, hover_data=index);
    # bokeh.plotting.output_file(html_path)
    # bokeh.plotting.save(p)
    # print(f'Saved plot HTML to: {html_path}')

def plot_umap_each_target(index_file, n_unit_threshold=1, cluster=False, min_samples=5, min_cluster_size=5, cmap="Spectral", width=800, height=800, show_legend=True):

    dir = os.path.dirname(index_file)
    index = pd.read_csv(index_file, sep='\t')
    targets = index["target"]
    unique_target = np.unique(targets)
    for target in unique_target:
        png_path = os.path.join(dir, target + ".png")
        subindex = index[index["target"] == target]
        subindex = remove_row_by_unit_occurance(subindex, n=n_unit_threshold)
        points = subindex[["umap1", "umap2"]].to_numpy()
        
        print(f'\n> Drawing PNG for {target}...') 
        ax = plot_points(points, labels=subindex['unit'], markers=subindex['source'], cmap=cmap, width=width, height=height)
        ax.figure.savefig(png_path, bbox_inches='tight')
        print(f'Saved PNG to: {png_path}')

        if cluster:
            print(f'\n> Clustering {target}...')
            labels, clustered = fit_hdbscan(points, min_samples, min_cluster_size)
            png_path = os.path.join(dir, target + "_cluster.png")
            plot_hdbscan(points, labels, clustered, png_path, cmap)


if __name__ == "__main__":
    index_file = ".\\..\\..\\test_umap\\umap_result\\index.tsv"
    plot_umap_each_target(index_file, n_unit_threshold=15, cluster=True, min_samples=5, min_cluster_size=5, cmap="rainbow")