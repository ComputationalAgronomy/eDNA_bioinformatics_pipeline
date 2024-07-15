import hdbscan
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

def fit_hdbscan(
        points: np.ndarray,
        min_samples: int,
        min_cluster_size: int,
        cluster_selection_epsilon: float,
        alpha: float
    ) -> tuple[np.ndarray, np.ndarray, int, float]:
    """
    Fit HDBSCAN clustering to the given points and return the cluster labels, whether each point is clustered, the number of clusters, and the percentage of clustered points.

    :param points: Input data points (from UMAP embeddings) for HDBSCAN clustering.
    :param min_samples: Provide a measure of how conservative want clustering to be for HDBSCAN.
    :param min_cluster_size: Minimum number of samples required to consider as a cluster for HDBSCAN.
    :param cluster_selection_epsilon: The distance threshold that clusters below the given value are not split up any further.
    :param alpha: Determines how conservative HDBSCAN will try to cluster points together. Higher values will make HDBSCAN more conservative.
    """
    labels = hdbscan.HDBSCAN(
        min_samples=min_samples,
        min_cluster_size=min_cluster_size,
        cluster_selection_epsilon=cluster_selection_epsilon,
        alpha=alpha,
        allow_single_cluster=True
    ).fit_predict(points)
    clustered = (labels >= 0)
    n_c = max(labels) + 1
    noise = sum(1 for i in labels if i < 0)
    c_p = (1 - noise / len(labels)) * 100
    return labels, clustered, n_c, round(c_p, 2)

def plot_hdbscan(
        points: np.ndarray,
        labels: np.ndarray,
        clustered: np.ndarray,
        png_path: str,
        cmap: str,
        background: str = "white",
        width: int = 800,
        height:int = 800
    ) -> None:
    """
    Plot the HDBSCAN clustering results using matplotlib.
    """
    dpi = plt.rcParams["figure.dpi"]
    fig = plt.figure(figsize=(width / dpi, height / dpi))
    ax = fig.add_subplot(111)
    ax.set_facecolor(background)
    
    point_size = 300.0 / np.sqrt(points.shape[0])

    ax.scatter(
        points[~clustered, 0],
        points[~clustered, 1],
        color=(0.5, 0.5, 0.5),
        s=point_size,
        alpha=0.5
    )
    ax.scatter(
        points[clustered, 0],
        points[clustered, 1],
        c=labels[clustered],
        s=point_size,
        cmap=cmap
    )
    ax.tick_params(bottom=False, left=False, labelbottom=False, labelleft=False)
    ax.figure.savefig(png_path)

def run_hdbscan(
        index: pd.DataFrame,
        min_samples: int,
        min_cluster_size: int,
        cluster_selection_epsilon: float,
        alpha: float,
        png_path: str,
        cmap: str
    ) -> tuple[str, str, str]:
    """
    run HDBSCAN clustering and plot the results.
    """
    points = index[["umap1", "umap2"]].to_numpy()
    labels, clustered, numb_clus, clus_perc = fit_hdbscan(points, min_samples, min_cluster_size, cluster_selection_epsilon, alpha)
    numb_unit = len(index["unit"].unique())
    plot_hdbscan(points, labels, clustered, png_path, cmap)
    print(f'Saved PNG to: {png_path}')
    return str(numb_unit), str(numb_clus), str(clus_perc)

def write_cluster_report(cluster_report: list[list[str]], prefix: str, save_dir: str) -> None:
    """
    Write the cluster report to a TSV file.
    Columns: name, unit_counts, cluster_counts, clustered_ratio
    """
    cluster_report_path = os.path.join(save_dir, f"{prefix}_cluster_report.tsv")
    cluster_report = pd.DataFrame(
        cluster_report,
        columns=["name", "unit_counts", "cluster_counts", "clustered_ratio"]
    )
    cluster_report.to_csv(cluster_report_path, sep='\t', index=False)
    print(f'Saved cluster report to: {cluster_report_path}')

def run_hdbscan_by_category(
        index: pd.DataFrame,
        min_samples: int,
        min_cluster_size: int,
        cluster_selection_epsilon: float,
        alpha: float,
        category: str,
        prefix: str,
        save_dir: str,
        cmap: str,
    ) -> None:
    """
    Run HDBSCAN clustering for a specified category and plot the results, grouped by the specified category.
    Writes a cluster report containing the number of units, clusters, and clustered ratio for each group.

    :param index: Input index DataFrame.
    :param min_samples: Provide a measure of how conservative want clustering to be for HDBSCAN.
    :param min_cluster_size: Minimum number of samples required to consider as a cluster for HDBSCAN.
    :param cluster_selection_epsilon: The distance threshold that clusters below the given value are not split up any further.
    :param alpha: Determines how conservative HDBSCAN will try to cluster points together. Higher values will make HDBSCAN more conservative.
    :param category: Column name to group the units by, restricted to 'unit', 'target', or 'all'.
    :param prefix: Prefix for the output file names.
    :param save_dir: Directory to save the output files.
    :param cmap: Color map for the plot.
    """
    if category == 'all':
        print("> Clustering for all units...")
        png_path = os.path.join(save_dir, f"{prefix}_hdbscan.png")

        numb_unit, numb_clus, clus_perc = run_hdbscan(
            index=index,
            min_samples=min_samples,
            min_cluster_size=min_cluster_size,
            cluster_selection_epsilon=cluster_selection_epsilon,
            alpha=alpha,
            png_path=png_path,
            cmap=cmap
        )

        cluster_report = ["all_targets", numb_unit, numb_clus, clus_perc]
        write_cluster_report(cluster_report, prefix, save_dir)

        return

    cluster_report = []

    unique_values = np.unique(index[category])
    for value in unique_values:
        print(f"> Clustering for {category}: {value}...")
        png_path = os.path.join(save_dir, f"{prefix}_{value}_hdbscan.png")

        subindex = index[index[category] == value]

        numb_unit, numb_clus, clus_perc = run_hdbscan(
            index=subindex,
            min_samples=min_samples,
            min_cluster_size=min_cluster_size,
            cluster_selection_epsilon=cluster_selection_epsilon,
            alpha=alpha,
            png_path=png_path,
            cmap=cmap
        )
        cluster_report.append([value, numb_unit, numb_clus, clus_perc])

    write_cluster_report(cluster_report, prefix, save_dir)    

    return