import argparse

def parse_args():
    level_choices = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    parser = argparse.ArgumentParser(prog="ednapc", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
    subparsers = parser.add_subparsers(help="sub-command help", metavar="SUBCOMMAND", required=True)

    # relative abundance barchart
    bc_parser = subparsers.add_parser('-barchart', help="Plot a barchart to visualize the relative abundance of a level across samples.")
    bc_parser.add_argument('-l', '--level',  choices=level_choices, type=str, required=True, help="The level to plot the barchart for.")
    bc_parser.add_argument('-s', '--save-dir', default='.', type=str, help="The directory to save the barchart HTML file.")
    bc_parser.add_argument('-p', '--prefix', default='barchart', type=str, help="The prefix for the barchart HTML file.")
    bc_parser.add_argument('-sil', '--sample-id-list', nargs='+', default=[], type=str, help="A list of sample IDs to plot. Default is None (plot all samples).")

    # maximum likelihood tree
    ml_parser = subparsers.add_parser('-tree', help="Reconstruct a phylogenetic tree for a list of targets using IQTREE.")
    ml_parser.add_argument('-tl', '--target-list', nargs='+', type=str, required=True, help="A list of targets to be plotted, e.g. FamilyA FamilyB FamilyC.")
    ml_parser.add_argument('-tlv', '--target-levl', choices=level_choices, type=str, required=True, help="The taxonomic level of the targets.")
    ml_parser.add_argument('-ulv', '--unit-levl', default='species', choices=level_choices,type=str, help="The taxonomic level of the units.")
    ml_parser.add_argument('-s', '--save-dir', default='.', type=str, help="The directory to save the output files.")
    ml_parser.add_argument('-p', '--prefix', default='mltree', type=str, help="The prefix for the output file names.")
    ml_parser.add_argument('-m', '--model', default=None, type=str, help="The model to specify for tree inference. If not specified, it will use the best-fit model found.")
    ml_parser.add_argument('-b', '--bootstrap', default=None, type=int, help="The number of bootstrap replicates. Default is None.")
    ml_parser.add_argument('-t', '--threads', default=None, type=int, help="The number of threads to use. If not specified, it will automatically determine the best number of cores given the current data and computer.")
    ml_parser.add_argument('-sil', '--sample-id-list', nargs='+', default=[], type=str, help="A list of sample IDs to plot. Default is None (plot all samples).")
    
    # UMAP
    umap_parser = subparsers.add_parser('-umap', help="Generate and plot UMAP visualizations for a list of targets based on their sequence data.")
    umap_parser.add_argument('-tl', '--target-list', nargs='+', type=str, required=True, help="A list of targets to be plotted, e.g. FamilyA FamilyB FamilyC.")
    umap_parser.add_argument('-tlv', '--target-levl', choices=level_choices, type=str, required=True, help="The taxonomic level of the targets.")
    umap_parser.add_argument('-ulv', '--unit-levl', default='species', choices=level_choices, type=str, help="The taxonomic level of the units.")
    umap_parser.add_argument('-nn', '--n-neighbors', default=15, type=int, help="The number of neighbors to consider for UMAP.")
    umap_parser.add_argument('-md', '--min-dist', default=0.1, type=float, help="The minimum distance parameter for UMAP.")
    umap_parser.add_argument('-rs', '--random-state', default=42, type=int, help="The random seed for reproducibility.")
    umap_parser.add_argument('-cd', '--calc-dist', action='store_true', help="If True, calculates a distance matrix for UMAP. Otherwise, transforms sequences into a one-hot encoded matrix.")
    umap_parser.add_argument('-s', '--save-dir', default='.', type=str, help="The directory to save the output files.")
    umap_parser.add_argument('--plot-all', action='store_true', help="If specified, plot all targets in one figure.")
    umap_parser.add_argument('--plot-target', action='store_true', help="If specified, plot each target in separate figure.")
    umap_parser.add_argument('--plot-unit', action='store_true', help="If specified, plot each unit in separate figure.")
    umap_parser.add_argument('-ip', '--index-path', default=None, type=str, help="Path to a pre-created index file. If provided, UMAP embedding calculation will be skipped.")
    umap_parser.add_argument('-ds', '--dereplicate-sequence', action='store_true', help="If specified, use unique sequences as input data for UMAP.")
    umap_parser.add_argument('--cmap', default='rainbow', type=str, help="The colormap for the plots.")
    umap_parser.add_argument('-nl', '--no-legend', action='store_false', help="If specified, don't show legend in the plots.")
    umap_parser.add_argument('-sil', '--sample-id-list', nargs='+', default=[], type=str, help="A list of sample IDs to plot. Default is None (plot all samples).")

    # cluster UMAP
    cu_parser = subparsers.add_parser('-cluster', help="Cluster UMAP embeddings and generate plots for a list of targets based on a pre-created umap index file.")
    cu_parser.add_argument('-ip', '--index-path', type=str, required=True, help="Path to the pre-created umap index file.")
    cu_parser.add_argument('-n', '--n-unit-threshold', type=int, required=True, help="Filter out units that occur less than n times in the index DataFrame. The value should be set equal to UMAP n_neighbors.")
    cu_parser.add_argument('-ms', '--min-samples', default=5, type=int, help="Provide a measure of how conservative want clustering to be for HDBSCAN.")
    cu_parser.add_argument('-mcs', '--min-cluster-size', default=5, type=int, help="Minimum number of samples required to consider as a cluster for HDBSCAN.")
    cu_parser.add_argument('-cse', '--cluster-selection-epsilon', default=1.0, type=float, help="The distance threshold that clusters below the given value are not split up any further.")
    cu_parser.add_argument('-a', '--alpha', default=1.0, type=float, help="Determines how conservative HDBSCAN will try to cluster points together. Higher values will make HDBSCAN more conservative.")
    cu_parser.add_argument('-s', '--save-dir', default='.', type=str, help="The directory to save the output files.")
    umap_parser.add_argument('--plot-all', action='store_true', help="If specified, plot all targets in one figure.")
    umap_parser.add_argument('--plot-target', action='store_true', help="If specified, plot each target in separate figure.")
    umap_parser.add_argument('--plot-unit', action='store_true', help="If specified, plot each unit in separate figure.")
    cu_parser.add_argument('--cmap', default='rainbow', type=str, help="The colormap for the plots.")

    # nexus file
    nex_parser = subparsers.add_parser('-nex', help="Generate NEXUS file for a given species. The NEXUS file contains the aligned sequences, their corresponding labels, and the frequency of each label for each unique sequence.")
    nex_parser.add_argument('-ip', '--index-path', type=str, required=True, help="Path to the index file containing the unit information.")
    nex_parser.add_argument('-sn', '--species-name', type=str, required=True, help="Name of the species for which the NEXUS file will be generated.")
    nex_parser.add_argument('-lt', '--label-type', choices=['hdbscan', 'site'], type=str, required=True, help="Type of labels to use. 'hdbscan' label type uses HDBSCAN clustering results for labeling and calculates the frequency of each cluster label. 'site' label type uses site information (e.g., 'taoyuan' or 'keelung') for labeling and calculates the frequency of each site label.")
    nex_parser.add_argument('-s', '--save-dir', default='.', type=str, help="Directory where the NEXUS file will be saved.")
    nex_parser.add_argument('sil', '--sample-id-list', nargs='+', default=[], type=str, help="A list of sample IDs to include in the NEXUS file. The list should be same as that specified by the index file. If empty, all samples will be included.")

    return parser.parse_args()

def main(args):
    pass

if __name__ == "__main__":
    args = parse_args()
    main(args)