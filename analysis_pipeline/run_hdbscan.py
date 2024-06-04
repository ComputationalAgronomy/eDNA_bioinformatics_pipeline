import hdbscan
import matplotlib.pyplot as plt
import numpy as np
from umap.plot import _themes

def fit_hdbscan(points, min_samples=5, min_cluster_size=5):
    labels = hdbscan.HDBSCAN(
        min_samples=min_samples,
        min_cluster_size=min_cluster_size
        ).fit_predict(points)
    clustered = (labels >= 0)
    num_c = max(labels) + 1
    noise = sum(1 for i in labels if i < 0)
    print(f'Number of clusters: {num_c}')
    print(f'clustered percentage: {100 - noise / len(labels) * 100:.2f}%')
    return labels, clustered

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

