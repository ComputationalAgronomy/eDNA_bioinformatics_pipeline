from analysis_function import align_fasta, create_dir
from color_umap import plot_points
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
            label = record.description.split("-")[0]
            seq_id = record.description
            index_list.append([index, label, seq_id])
            record.id = str(i)
            record.description = ''
            record.name = record.id
            SeqIO.write(record, file, 'fasta')
    
    index = pd.DataFrame(index_list, columns=["index", "label", "seq_id"])
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

def run_umap(seq_file, umap_dir, target_name, random_state=1, neighbors=15, min_dist=0.1, theme='fire', width=800, height=800):
    target_dir = os.path.join(umap_dir, target_name)
    create_dir(target_dir)
    index_fasta_file = os.path.join(target_dir, "input.fa")
    aln_file = os.path.join(target_dir, f"input.aln")
    dist_file = os.path.join(target_dir, "distance.txt")
    index_path = os.path.join(target_dir, "index.tsv")
    png_path = os.path.join(target_dir, "umap.png") 
    html_path = os.path.join(target_dir, "umap.html")

    index = fasta2index(seq_file=seq_file, fasta_file=index_fasta_file)
    calc_dist(seq_file=index_fasta_file, aln_file=aln_file, dist_file=dist_file)
    dist_matrix = load_sparse_dist_matrix(dist_file)
    reducer, embedding = fit_umap(dist_matrix, random_state=random_state, neighbors=neighbors, min_dist=min_dist)
    index['umap1'] = embedding[:,0]
    index['umap2'] = embedding[:,1]

    index.to_csv(index_path, sep='\t', index=False)
    print(f'Saved index TSV to: {index_path}')

    print('\n> Drawing PNG...')
    ax = plot_points(reducer, labels=index['label'], theme=theme, width=width, height=height)
    ax.figure.savefig(png_path, bbox_inches='tight')
    print(f'Saved PNG to: {png_path}')
    
    print('\n> Drawing interactive plot...')
    p = umap.plot.interactive(reducer, labels=index['label'], theme=theme, width=width, height=height, hover_data=index);
    bokeh.plotting.output_file(html_path)
    bokeh.plotting.save(p)
    print(f'Saved plot HTML to: {html_path}')
    
    print(f'\nDone. Saved to: {target_dir}')
    return reducer, index