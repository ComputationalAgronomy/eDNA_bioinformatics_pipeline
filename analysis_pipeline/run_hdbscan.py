import hdbscan
import matplotlib.pyplot as plt
import numpy as np

def fit_hdbscan(points, min_samples=5, min_cluster_size=5):
    labels = hdbscan.HDBSCAN(
        min_samples=min_samples,
        min_cluster_size=min_cluster_size,
        cluster_selection_epsilon=0.5,
        alpha=1.0,
        allow_single_cluster=True
        ).fit_predict(points)
    clustered = (labels >= 0)
    n_c = max(labels) + 1
    noise = sum(1 for i in labels if i < 0)
    c_p = (1 - noise / len(labels)) * 100
    return labels, clustered, n_c, round(c_p, 2)

def plot_hdbscan(points, labels, clustered, png_path, cmap="Spectral", background="white", width=800, height=800):
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
        cmap=cmap)
    ax.figure.savefig(png_path)

def run_hdbscan(index, min_samples=5, min_cluster_size=5, png_path="./clustered.png", cmap="Spectral"):
    print(f'\n> Clustering...')
    points = index[["umap1", "umap2"]].to_numpy()
    labels, clustered, numb_clus, clus_perc = fit_hdbscan(points, min_samples, min_cluster_size)
    numb_unit = len(index["unit"].unique())
    plot_hdbscan(points, labels, clustered, png_path, cmap)
    print(f'Saved PNG to: {png_path}')
    return numb_unit, numb_clus, clus_perc
