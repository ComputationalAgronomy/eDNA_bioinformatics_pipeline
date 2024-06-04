import hdbscan
import pandas as pd
import matplotlib.pyplot as plt
from run_umap import remove_row_by_unit_occurance

def fit_hdbscan(points, min_samples=5, min_cluster_size=5):
    labels = hdbscan.HDBSCAN(
        min_samples=min_samples,
        min_cluster_size=min_cluster_size
        ).fit_predict(points)
    clustered = (labels >= 0)
    plt.scatter(
        points[~clustered, 0],
        points[~clustered, 1],
        color=(0.5, 0.5, 0.5),
        s=0.1,
        alpha=0.5
        )
    plt.scatter(
        points[clustered, 0],
        points[clustered, 1],
        c=labels[clustered],
        s=0.1,
        cmap='Spectral');   
    plt.show()

if __name__ == "__main__":
    index_file = ".\\..\\..\\test_umap\\umap_result\\index.tsv"
    n_unit_threshold = 5
    index = pd.read_csv(index_file, sep='\t')
    index = remove_row_by_unit_occurance(index, n=n_unit_threshold)
    points = index[["umap1", "umap2"]].to_numpy()
    fit_hdbscan(points)
