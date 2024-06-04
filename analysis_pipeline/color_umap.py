import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.cm
from matplotlib.patches import Patch
from umap.plot import _themes, _datashade_points
import matplotlib.pyplot as plt
import os

def _matplotlib_points(
    points,
    ax=None,
    labels=None,
    markers=None,
    values=None,
    cmap="Blues",
    color_key=None,
    color_key_cmap="Spectral",
    background="white",
    width=800,
    height=800,
    show_legend=True,
    alpha=None,
    symbol_map = ["o", "D", "*", "s", "h", "8", "X", "p"]
):
    """Use matplotlib to plot points"""
    point_size = 300.0 / np.sqrt(points.shape[0])

    legend_elements = None

    if ax is None:
        dpi = plt.rcParams["figure.dpi"]
        fig = plt.figure(figsize=(width / dpi, height / dpi))
        ax = fig.add_subplot(111)

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
            color_key = plt.get_cmap(color_key_cmap)(np.linspace(0, 1, num_labels))
            legend_elements = [
                Patch(facecolor=color_key[i], label=unique_labels[i])
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

def plot_points(points, labels=None, markers=None, values=None, theme=None, cmap="Blues", color_key=None, color_key_cmap="Spectral", background="white", width=800, height=800, show_legend=True):
    if theme is not None:
        cmap = _themes[theme]["cmap"]
        color_key_cmap = _themes[theme]["color_key_cmap"]
        if background is None:
            background = _themes[theme]["background"]

    dpi = plt.rcParams["figure.dpi"]
    fig = plt.figure(figsize=(width / dpi, height / dpi))
    ax = fig.add_subplot(111)

    if points.shape[0] <= width * height // 10:
        ax = _matplotlib_points(points, ax, labels, markers, values, cmap, color_key, color_key_cmap, background, width, height, show_legend)
    else:
        ax = _datashade_points(points, ax, labels, values, cmap, color_key, color_key_cmap, background, width, height, show_legend)

    ax.set(xticks=[], yticks=[])
 
    return ax

def plot_by_index(index_file, values=None, theme=None, cmap="Blues", color_key=None, color_key_cmap="Spectral", background="white", width=800, height=800, show_legend=True):
    if theme is not None:
        cmap = _themes[theme]["cmap"]
        color_key_cmap = _themes[theme]["color_key_cmap"]
        if background is None:
            background = _themes[theme]["background"]

    dir = os.path.dirname(index_file)
    index = pd.read_csv(index_file, sep='\t')
    familys = index["family"]
    unique_family = np.unique(familys)
    for family in unique_family:
        print(family)
        png_path = os.path.join(dir, family + ".png")
        subindex = index[index["family"] == family]
        points = subindex[["umap1", "umap2"]].to_numpy()
        dpi = plt.rcParams["figure.dpi"]
        fig = plt.figure(figsize=(width / dpi, height / dpi))
        ax = fig.add_subplot(111)
        labels = subindex["label"]
        markers = subindex["source"]

        if points.shape[0] <= width * height // 10:
            ax = _matplotlib_points(points, ax, labels, markers, values, cmap, color_key, color_key_cmap, background, width, height, show_legend)
        else:
            ax = _datashade_points(points, ax, labels, values, cmap, color_key, color_key_cmap, background, width, height, show_legend)

        ax.set(xticks=[], yticks=[])
        ax.figure.savefig(png_path, bbox_inches='tight')

if __name__ == "__main__":
    index_file = ".\\test_umap\\species_index_larger_15.tsv"
    plot_by_index(index_file)