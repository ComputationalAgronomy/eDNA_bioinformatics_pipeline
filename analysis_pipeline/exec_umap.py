import os
import pandas as pd
from utils_umap import (
    fasta2index,
    compute_embedding_with_distmx,
    compute_embedding_with_onehotmx,
    update_index,
    plot_points
)
from utils_sequence import align_fasta

def write_umap_file(
        seq_path: str,
        save_dir: str,
        unit2target: dict[str, str]=None,
        random_state: int = 42,
        calc_dist: bool = True,
        neighbors: int = 15,
        min_dist: float = 0.1
    ) -> pd.DataFrame:
    """
    Run the UMAP pipeline and save the index TSV file.

    :param seq_path: Path to the input aligned FASTA file.
    :param save_dir: Path to the output directory.
    :param unit2target: Dictionary mapping unit labels to target labels, used to create the 'target' column for the index.. Default is None.
    :param random_state: Random state for umap. Defaults is 42.
    :param calc_dist: Whether to calculate the distance matrix as the input of umap. Default is True
    :param neighbors: Number of neighbors for umap. Defaults is 15.
    :param min_dist: Minimum distance for umap. Defaults is 0.1.
    :return: Index DataFrame.
    """
    os.makedirs(save_dir, exist_ok=True)

    index_fasta_path = os.path.join(save_dir, "input.fa")
    aln_path = os.path.join(save_dir, "input.aln")
    index_path = os.path.join(save_dir, "umap_index.tsv")

    index = fasta2index(seq_path=seq_path, index_fasta_path=index_fasta_path)

    align_fasta(seq_path=seq_path, aln_path=aln_path)

    if calc_dist:
        embedding = compute_embedding_with_distmx(aln_path, save_dir, random_state, neighbors, min_dist)
    else:
        embedding = compute_embedding_with_onehotmx(aln_path, random_state, neighbors, min_dist)

    index = update_index(index, unit2target, embedding)
    index.to_csv(index_path, sep='\t', index=False)
    print(f'Saved index TSV to: {index_path}')

    return index

def plot_umap(
        index,
        png_path='umap.png', 
        cmap="Spectral",
        width=800,
        height=800,
        show_legend=True
    ) -> None:
    """
    Plot the UMAP embedding and save the plot as a PNG file.

    :param index: Index DataFrame.
    :param png_path: Path to the output PNG file. Default is 'umap.png'.
    :param cmap: Color map for the plot. Default is 'Spectral'.
    :param width: Width of the plot in pixels. Default is 800.
    :param height: Height of the plot in pixels. Default is 800.
    :param show_legend: Whether to show the legend. Default is True.
    """
    points = index[["umap1", "umap2"]].to_numpy()
    print('\n> Drawing PNG...')
    ax = plot_points(
        points,
        labels=index['unit'],
        markers=index['source'],
        cmap=cmap,
        width=width,
        height=height,
        show_legend=show_legend
    )
    ax.figure.savefig(png_path, bbox_inches='tight')
    print(f'Saved PNG to: {png_path}')

    # print('\n> Drawing interactive plot...')
    # p = umap.plot.interactive(reducer, labels=index['label'], theme=theme, width=width, height=height, hover_data=index);
    # bokeh.plotting.output_file(html_path)
    # bokeh.plotting.save(p)
    # print(f'Saved plot HTML to: {html_path}')