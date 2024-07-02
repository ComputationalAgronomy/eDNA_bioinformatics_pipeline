from analysis_function import align_fasta, create_dir
from color_umap import plot_points
import os
import numpy as np
import pandas as pd
from scipy import sparse
import subprocess
from Bio import SeqIO
import umap

def fasta2index(seq_file, fasta_file):
    index_list = []
    with open(fasta_file, 'w')as file:
        records = SeqIO.parse(seq_file, 'fasta')
        for i, record in enumerate(records):
            index = str(i)
            unit = record.description.rsplit("-", 1)[0]
            seq_id = record.description
            index_list.append([index, seq_id, unit])
            record.id = str(i)
            record.description = ''
            record.name = record.id
            SeqIO.write(record, file, 'fasta')
    
    index = pd.DataFrame(index_list, columns=["index", "seq_id", "unit"])
    return index

def calc_dist(seq_file, dist_file, maxdist=1.0, termdist=1.0, threads=12):
    cmd = f"usearch -calc_distmx {seq_file} -tabbedout {dist_file} -maxdist {str(maxdist)} -termdist {str(termdist)}"
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

def sequence2number(sequence):
    base_map = {'A':[1, 0, 0, 0], 
                'C':[0, 1, 0, 0], 
                'G':[0, 0, 1, 0], 
                'T':[0, 0, 0, 1]}
    number = []
    for base in sequence:
        number.extend(base_map.get(base, [0, 0, 0, 0]))
    return number

def create_number_matrix(seq_file):
    number_matrix = []
    with open(seq_file, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            number_matrix.append(sequence2number(record.seq))
    return np.array(number_matrix)

def fit_umap(matrix, random_state=None, neighbors=15, min_dist=0.1, spread=1.0, precomputed=True):
    print(f'> Creating UMAP embedding with {neighbors} neighbors...')
    reducer = umap.UMAP(n_neighbors=neighbors, random_state=random_state, min_dist=min_dist, spread=spread, metric="precomputed" if precomputed else "euclidean")
    embedding = reducer.fit_transform(matrix)
    return reducer, embedding

def get_index_source_label(seq_id):
    s = []
    for id in seq_id:
        if "taoyuan" in id:
            s.append("taoyuan")
        elif "keelung" in id:
            s.append("keelung")
        else:
            s.append("unknown")
    return s

def get_index_target_label(index, target2units):
    t = []
    for unit in index['unit']:
        for target, units in target2units.items():
            if unit in units:
                t.append(target)
    return t

def write_umap_file(seq_file, save_dir, target2units=None, random_state=42, calc_dist=True, neighbors=15, min_dist=0.1):
    create_dir(save_dir)
    index_fasta_file = os.path.join(save_dir, "input.fa")
    aln_file = os.path.join(save_dir, "input.aln")
    index_path = os.path.join(save_dir, "index.tsv")

    index = fasta2index(seq_file=seq_file, fasta_file=index_fasta_file)
    if target2units is not None:
        index['target'] = get_index_target_label(index, target2units)
    index['source'] = get_index_source_label(index['seq_id'])

    align_fasta(seq_file=index_fasta_file, aln_file=aln_file)

    if calc_dist:
        dist_file = os.path.join(save_dir, "distance.txt")
        calc_dist(seq_file=aln_file, dist_file=dist_file)
        dist_matrix = load_sparse_dist_matrix(dist_file)
        reducer, embedding = fit_umap(dist_matrix, random_state=random_state, neighbors=neighbors, min_dist=min_dist, precomputed=True)
    else:
        number_matrix = create_number_matrix(aln_file)
        reducer, embedding = fit_umap(number_matrix, random_state=random_state, neighbors=neighbors, min_dist=min_dist, precomputed=False)
    index['umap1'] = embedding[:,0]
    index['umap2'] = embedding[:,1]

    index.to_csv(index_path, sep='\t', index=False)
    print(f'Saved index TSV to: {index_path}')
    return index

def plot_umap(index, png_path='./umap.png', cmap="Spectral", width=800, height=800, show_legend=True):
    points = index[["umap1", "umap2"]].to_numpy()
    print('\n> Drawing PNG...')
    ax = plot_points(points, labels=index['unit'], markers=index['source'], cmap=cmap, width=width, height=height, show_legend=show_legend)
    ax.figure.savefig(png_path, bbox_inches='tight')
    print(f'Saved PNG to: {png_path}')

    # print('\n> Drawing interactive plot...')
    # p = umap.plot.interactive(reducer, labels=index['label'], theme=theme, width=width, height=height, hover_data=index);
    # bokeh.plotting.output_file(html_path)
    # bokeh.plotting.save(p)
    # print(f'Saved plot HTML to: {html_path}')
