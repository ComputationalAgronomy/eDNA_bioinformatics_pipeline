from Bio import SeqIO
import matplotlib.cm
import matplotlib.colors
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy import sparse
import subprocess
import tempfile
import umap
from umap.plot import _datashade_points, _themes

from edna_processor.utils.base_logger import logger
from edna_processor.utils.utils_sequence import write_fasta, align_fasta

def fasta2index(seq_path: str, index_fasta_path: str) -> pd.DataFrame:
    """
    Read in a FASTA file and output an index FASTA file replacing the sequence IDs with indexes.

    :param seq_path: Path to the input FASTA file.
    :param index_fasta_path: Path to the output index FASTA file.
    :return: DataFrame containing the index, sequence ID, and unit name.
    """
    index_list = []

    with open(seq_path, 'r') as in_handle, open(index_fasta_path, 'w') as out_handle:
        for i, record in enumerate(SeqIO.parse(in_handle, 'fasta')):
            index = str(i)
            unit = record.description.rsplit("-", 1)[0]
            seq_id = record.description

            index_list.append([index, seq_id, unit])

            record.id = index
            record.description = ''
            record.name = index

            SeqIO.write(record, out_handle, 'fasta')

    index_df = pd.DataFrame(index_list, columns=["index", "seq_id", "unit"])

    return index_df

def get_index_source_label(seq_id: list[str]) -> list[str]:
    """
    Take a list of sequence IDs and return a list of source labels. 
    It identifies if the sequence ID contains the substring recorded in the 'sources' variable and assigns the corresponding label. 
    If neither substring is found, it assigns the label "unknown".

    :param seq_id: A list of sequence IDs.
    :return: A list of source labels corresponding to each sequence ID.
    """
    sources = ["taoyuan", "keelung"]

    source_labels = []
    for id in seq_id:
        label = "unknown"
        for source in sources:
            if source in id:
                label = source
                break
        source_labels.append(label)

    return source_labels

def get_index_target_label(unit_labels: list[str], unit2target: dict[str, str]) -> list[str]:
    """
    Take a list of unit labels and return a list of target labels.
    It uses the 'unit2target' dictionary to map unit labels to target labels. 
    If a unit label is not found in the dictionary, it assigns the label "unknown".

    :param unit_labels: A list of unit labels.
    :param unit2target: A dictionary mapping unit labels to target labels.
    :return: A list of target labels corresponding to each unit label.
    """
    target_labels = []
    for unit in unit_labels:
        if unit in unit2target:
            target_labels.append(unit2target[unit])
        else:
            target_labels.append("unknown")

    return target_labels

def calc_distmx(seq_path: str, dist_path: str, maxdist: float = 1.0, termdist: float = 1.0, threads: int = 12) -> None:
    """
    Calculate distance matrix using USEARCH.
    (USEARCH command reference: https://drive5.com/usearch/manual/cmd_calc_distmx.html)

    :param seq_path: Path to the input aligned FASTA file.
    :param dist_path: Path to the output distance matrix file.
    :param maxdist: The maximum distance to be written. Default is 1.0.
    :param termdist: The distance threshold for terminating the calculation. Default is 1.0.
    :param threads: Number of threads to use for the calculation. Default is 12.
    """
    print("> Calculating distance matrix...")

    cmd = [
        "usearch", "-calc_distmx", seq_path, "-tabbedout", dist_path,
        "-maxdist", str(maxdist), "-termdist", str(termdist)
    ]
    if threads:
        cmd.extend(["-threads", str(threads)])

    logger.info("> Running USEARCH command:", cmd)
    try:
        subprocess.run(cmd, check=True)
        logger.info(f"USEARCH finished. Output distance matrix file saved to: {dist_path}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred during the calculation of the distance matrix: {e}")

def load_sparse_dist_matrix(dist_path: str) -> np.ndarray:
    """
    Load a sparse distance matrix from a distance matrix file created by the 'calc_distmx' function.

    :param dist_path: Path to the input distance matrix file.
    :return: Sparse distance matrix as a NumPy array.
    """
    dist_matrix = pd.read_csv(dist_path, header=None, sep='\t')
    logger.info(f"Loading sparse {max(dist_matrix[0])+1} x {max(dist_matrix[0])+1} distance matrix from: {dist_path}")
    
    diagonal = dist_matrix[0] == dist_matrix[1]
    row = np.concatenate([dist_matrix[0], dist_matrix[1][~diagonal]])
    col = np.concatenate([dist_matrix[1], dist_matrix[0][~diagonal]])
    data = 1 - np.concatenate([dist_matrix[2], dist_matrix[2][~diagonal]])

    dist_matrix = sparse.csr_matrix((data, (row, col)), dtype=np.float32)

    return 1 - dist_matrix.toarray()

def sequence_to_one_hot(sequence: str) -> list[int]:
    """
    Convert a sequence with only ATCG bases to a one-hot encoded vector.

    :param sequence: DNA sequence.
    :return: One-hot encoded vector.
    """
    base_map = {'A':[1, 0, 0, 0], 
                'C':[0, 1, 0, 0], 
                'G':[0, 0, 1, 0], 
                'T':[0, 0, 0, 1]}

    one_hot_encoded = []
    for base in sequence:
        one_hot_encoded.extend(base_map.get(base, [0, 0, 0, 0]))

    return one_hot_encoded

def create_one_hot_matrix(seq_path: str) -> np.ndarray:
    """
    Read in a aligned FASTA file and output a one-hot encoded matrix.

    :param seq_path: Path to the input FASTA file.
    :return: One-hot encoded matrix as a NumPy array.
    """
    logger.info(f"Creating one-hot encoded matrix from: {seq_path}")

    one_hot_matrix = []
    with open(seq_path, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            one_hot_matrix.append(sequence_to_one_hot(record.seq))

    return np.array(one_hot_matrix)

def fit_umap(
        matrix: np.ndarray,
        neighbors: int,
        min_dist: float,
        random_state: int,
        precomputed: bool,
    ) -> tuple[umap.UMAP, np.ndarray]:
    """
    Fit UMAP and return the UMAP object and the embedding.

    :param matrix: distance matrix or one-hot encoded matrix
    :param neighbors: Number of neighbors for umap. 
    :param min_dist: Minimum distance for umap. 
    :param random_state: Random state for umap. 
    :param precomputed: Whether the elements of the matrix are distances or not. 
    :return: UMAP object and embedding
    """
    logger.info(f'Creating UMAP embedding with {neighbors} neighbors...')

    reducer = umap.UMAP(
        n_neighbors=neighbors,
        min_dist=min_dist,
        random_state=random_state,
        metric="precomputed" if precomputed else "euclidean"
    )

    embedding = reducer.fit_transform(matrix)

    return reducer, embedding

def compute_embedding_with_distmx(
        seq_path: str,
        save_distmx_dir: str,
        neighbors: int,
        min_dist: float,
        random_state: int,
    ) -> np.ndarray:
    """
    The steps for calculating distance matrix and fitting UMAP.
    """
    dist_path = os.path.join(save_distmx_dir, "distance.txt")
    calc_distmx(seq_path=seq_path, dist_path=dist_path)
    dist_matrix = load_sparse_dist_matrix(dist_path)
    reducer, embedding = fit_umap(
        matrix=dist_matrix,
        neighbors=neighbors,
        min_dist=min_dist,
        random_state=random_state,
        precomputed=True
    )
    return embedding

def compute_embedding_with_onehotmx(
        seq_path: str,
        neighbors: int,
        min_dist: float,
        random_state: int,
    ) -> np.ndarray:
    """
    The steps for fitting UMAP with one-hot encoded matrix.
    """
    number_matrix = create_one_hot_matrix(seq_path)
    reducer, embedding = fit_umap(
        number_matrix,
        neighbors=neighbors,
        min_dist=min_dist,
        random_state=random_state,
        precomputed=False
    )
    return embedding

def update_index(
        index: pd.DataFrame,
        unit2target: dict[str, str],
        embedding: np.ndarray
    ) -> pd.DataFrame:
    """
    The steps for updating the index DataFrame with source/target labels and UMAP cordinates.
    """
    if unit2target is not None:
        index['target'] = get_index_target_label(index["unit"], unit2target)
    index['source'] = get_index_source_label(list(index['seq_id']))
    index["umap1"] = embedding[:,0]
    index["umap2"] = embedding[:,1]

    return index

def filter_index_by_unit_occurrence(index: pd.DataFrame, n: int = 1) -> pd.DataFrame:
    """
    Filter out units that occur less than n times in the index DataFrame.
    The value of 'n' should be set equal to the value of the UMAP 'neighbor'.

    :param index: Index DataFrame.
    :param n: The threshold of minimum occurrence to keep a unit.
    :return: Filtered index DataFrame.
    """
    if n <= 1:
        return index
    counts = index["unit"].value_counts()
    units_to_remove = counts[counts < n].index
    filtered_index = index[~index["unit"].isin(units_to_remove)]
    logger.info(f"Units with less than {n} occurrences have been removed.")
    return filtered_index

def _matplotlib_points(
    points,
    ax=None,
    labels=None,
    markers=None,
    values=None,
    color_key=None,
    cmap="Spectral",
    background="white",
    width=800,
    height=800,
    show_legend=True,
    alpha=None,
    symbol_map = ["o", "D", "*", "s", "h", "8", "X", "p"]
):

    point_size = 300.0 / np.sqrt(points.shape[0])

    legend_elements = None

    # if ax is None:
    #     dpi = plt.rcParams["figure.dpi"]
    #     fig = plt.figure(figsize=(width / dpi, height / dpi))
    #     ax = fig.add_subplot(111)

    ax.set_facecolor(background)

    # Color by labels
    if labels is not None:
        if labels.shape[0] != points.shape[0]:
            raise ValueError(
                "Labels must have a label for "
                "each sample (size mismatch: {} {})".format(
                    labels.shape[0], points.shape[0]
                )
            )
        if color_key is None:
            unique_labels = np.unique(labels)
            num_labels = unique_labels.shape[0]
            color_key = plt.get_cmap(cmap)(np.linspace(0, 1, num_labels))
            legend_elements = [
                Patch(facecolor=color_key[i], label=k)
                for i, k in enumerate(unique_labels)
            ]

        if isinstance(color_key, dict):
            colors = pd.Series(labels).map(color_key)
            unique_labels = np.unique(labels)
            legend_elements = [
                Patch(facecolor=color_key[k], label=k) for k in unique_labels
            ]
        else:
            unique_labels = np.unique(labels)
            if len(color_key) < unique_labels.shape[0]:
                raise ValueError(
                    "Color key must have enough colors for the number of labels"
                )
            new_color_key = {
                k: matplotlib.colors.to_hex(color_key[i])
                for i, k in enumerate(unique_labels)
            }
            legend_elements = [
                Patch(facecolor=color_key[i], label=k)
                for i, k in enumerate(unique_labels)
            ]
            colors = pd.Series(labels).map(new_color_key)

        if markers is not None:
            m = []
            unique_markers = np.unique(markers)
            if len(unique_markers) > len(symbol_map):
                raise ValueError(
                    "Too many unique markers for the number of labels, please customize 'symbol_map'."
                )
            for marker in markers:
                for i, k in enumerate(unique_markers):
                    if marker == k:
                        m.append(symbol_map[i])
        colors = list(colors)
        for i in range(len(points[:, 0])):
            ax.scatter(points[i, 0], points[i, 1], s=point_size, c=colors[i], marker=m[i], alpha=alpha)
        # ax.scatter(points[:, 0], points[:, 1], s=point_size, c=colors, markers=m, alpha=alpha)


    # Color by values
    elif values is not None:
        if values.shape[0] != points.shape[0]:
            raise ValueError(
                "Values must have a value for "
                "each sample (size mismatch: {} {})".format(
                    values.shape[0], points.shape[0]
                )
            )
        ax.scatter(
            points[:, 0], points[:, 1], s=point_size, c=values, cmap=cmap, alpha=alpha
        )

    # No color (just pick the midpoint of the cmap)
    else:

        color = plt.get_cmap(cmap)(0.5)
        ax.scatter(points[:, 0], points[:, 1], s=point_size, c=color)

    if show_legend and legend_elements is not None:
        ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5))

    return ax

def plot_points(points, labels=None, markers=None, values=None, color_key=None, cmap="rainbow", show_legend=True, background="white", width=800, height=800):

    dpi = plt.rcParams["figure.dpi"]
    fig = plt.figure(figsize=(width / dpi, height / dpi))
    ax = fig.add_subplot(111)

    if points.shape[0] <= width * height // 10:
        ax = _matplotlib_points(points, ax, labels, markers, values, color_key, cmap, background, width, height, show_legend)
    else:
        ax = _datashade_points(points, ax, labels, values, color_key, cmap, background, width, height, show_legend)

    ax.set(xticks=[], yticks=[])

    return ax

def run_umap(
        units2fasta: dict[str, str],
        dereplicate_sequence: bool,
        save_dir: str,
        neighbors: int,
        min_dist: float,
        random_state: int,
        calc_dist: bool,
        unit2target: dict[str, str]
    ) -> pd.DataFrame:
    """
    Run the UMAP pipeline and save the index TSV file.
    Step:
        1. Write the units FASTA to a file, with the sequence title formatted as >{unit_name}-{sample_id}_{haplotype_id}.
        2. Convert the unit FASTA file to as index FASTA file, with sequence title formatted as >{index}, and an index TSV file.
        3. Align the sequences using Clustal Omega.
        4. Compute the UMAP embedding using either the distance matrix or one-hot encoding.
        5. Update the index TSV file with the UMAP coordinates and target labels (if provided).
        6. Save the index TSV file to the specified output directory.
    """
    os.makedirs(save_dir, exist_ok=True)

    temp_dir = tempfile.TemporaryDirectory()
    unit_fasta_path = os.path.join(temp_dir.name, 'umap.fa')
    index_fasta_path = os.path.join(temp_dir.name, "input.fa")
    aln_fasta_path = os.path.join(save_dir, "input.aln")
    index_path = os.path.join(save_dir, "umap_index.tsv")

    write_fasta(units2fasta_dict=units2fasta, save_path=unit_fasta_path, dereplicate=dereplicate_sequence)

    index = fasta2index(seq_path=unit_fasta_path, index_fasta_path=index_fasta_path)

    align_fasta(seq_path=index_fasta_path, aln_path=aln_fasta_path)

    temp_dir.cleanup()

    if calc_dist:
        embedding = compute_embedding_with_distmx(
            seq_path=aln_fasta_path,
            save_distmx_dir=save_dir,
            neighbors=neighbors,
            min_dist=min_dist,
            random_state=random_state
        )
    else:
        embedding = compute_embedding_with_onehotmx(
            seq_path=aln_fasta_path,
            neighbors=neighbors,
            min_dist=min_dist,
            random_state=random_state
        )

    index = update_index(index, unit2target, embedding)
    index.to_csv(index_path, sep='\t', index=False)
    logger.info(f'Saved index TSV to: {index_path}')

    return index

def plot_umap(
        index: pd.DataFrame,
        png_path: str = 'umap.png', 
        cmap: str = 'rainbow',
        show_legend: bool = True
    ) -> None:
    """
    Plot the UMAP embedding and save the plot as a PNG file.
    """
    points = index[["umap1", "umap2"]].to_numpy()
    ax = plot_points(
        points=points,
        labels=index['unit'],
        markers=index['source'],
        cmap=cmap,
        show_legend=show_legend
    )
    ax.figure.savefig(png_path, bbox_inches='tight')
    logger.info(f"Saved PNG to: {png_path}")

    # print('\n> Drawing interactive plot...')
    # p = umap.plot.interactive(reducer, labels=index['label'], theme=theme, width=width, height=height, hover_data=index);
    # bokeh.plotting.output_file(html_path)
    # bokeh.plotting.save(p)
    # print(f'Saved plot HTML to: {html_path}')

def plot_umap_by_category(
        index: pd.DataFrame,
        category: str,
        prefix: str,
        png_dir: str,
        cmap: str,
        show_legend: bool
    ) -> None:
    """
    Plot the UMAP embedding and save the plot as a PNG file, grouped by the specified category.
    
    :param index: The index DataFrame containing the UMAP coordinates and category labels.
    :param category: Column name to group the units by, restricted to 'unit', 'target', or 'all'.
    :param prefix: The prefix for the output file names.
    :param png_dir: The directory to save the PNG files in.
    :param cmap: The colormap to use for the plots.
    :param show_legend: Whether to show the legend in the plots.
    """
    if category == 'all':
        logger.info("Drawing PNG for all units...")
        png_path = os.path.join(png_dir, f"{prefix}_umap.png")
        plot_umap(index=index, png_path=png_path, cmap=cmap, show_legend=show_legend)
        return

    unique_values = np.unique(index[category])
    for value in unique_values:
        logger.info(f"Drawing PNG for {category} {value}...")
        png_path = os.path.join(png_dir, f"{prefix}_{value}_umap.png")
        subindex = index[index[category] == value]
        plot_umap(index=subindex, png_path=png_path, cmap=cmap, show_legend=show_legend)