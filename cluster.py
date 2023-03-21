#!/usr/bin/env python3

import os
import argparse
from collections import defaultdict
import numba
import pandas as pd
import numpy as np
import fastcluster
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster, to_tree, dendrogram
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator, FuncFormatter


@numba.jit(nopython=True, nogil=False)
def squareform_matrix(X):
    """
    :param X: ndarray
    :return: ndarray
    """
    m, n = X.shape
    dm = np.empty((m, m), dtype=np.double)
    for i in range(0, m):
        for j in range(i, m):
            dm[i, j] = dm[j, i] = np.sum(X[i] != X[j])
    return dm


class PairwiseDistanceMatrix:
    def __init__(self, profile):
        self.transform(profile)

    def __call__(self):
        data = squareform_matrix(self.profile.values)
        return pd.DataFrame(data, index=self.profile.index, columns=self.profile.index)

    @staticmethod
    def _integer_encoding(s):
        mapper = defaultdict(int)
        unique_codes = set(s.dropna())
        mapper.update({j: i for i, j in enumerate(unique_codes, 1)})
        s = s.map(mapper)
        return s

    def transform(self, profile):
        self.profile = profile.apply(self._integer_encoding, axis=1).T


class HierarchicalCluster:
    def __init__(self, distmatrix):
        self.cluster = fastcluster.single(squareform(distmatrix))


class Figure:
    plt.style.use("fast")
    rcParams["lines.linewidth"] = 0.5
    rcParams["svg.fonttype"] = "none"

    def __init__(self, width, height):
        self.fig, self.ax = plt.subplots(1, 1, figsize=(width, height))
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['bottom'].set_visible(False)
        self.ax.spines['left'].set_visible(False)
        self.ax.spines['right'].set_visible(False)
        self.ax.grid(False)
        self.ax.patch.set_facecolor('none')
        plt.close()

    def annotate(self, text, position, fontsize=8):
        self.ax.annotate(text, position, xytext=(-2, 8), textcoords='offset points', va='top', ha='right',
                         fontsize=fontsize)

    def savefig(self, file, dpi=100):
        self.fig.savefig(file, dpi=dpi, bbox_inches='tight', pad_inches=.5)


class Dendrogram:
    def __init__(self, profile):
        distmatrix = PairwiseDistanceMatrix(profile)()
        self.hc = HierarchicalCluster(distmatrix)
        self.figure = Figure(10, distmatrix.shape[1] * 0.3)
        self._labels = distmatrix.columns

    def savefig(self, outfile, dpi=100):
        self.figure.savefig(file=outfile, dpi=dpi)

    def to_newick(self, outfile):
        tree = to_tree(self.hc.cluster, False)
        newick = make_newick(tree, "", tree.dist, self._labels)
        with open(outfile, 'w') as handle:
            handle.write(newick)

    def __call__(self, no_labels=False, show_node_info=False, labels_color=None, xlim=None):
        tree = dendrogram(
            self.hc.cluster,
            ax=self.figure.ax,
            labels=self._labels,
            orientation="left",
            leaf_font_size=12,
            above_threshold_color="#000000",
            color_threshold=0,
            no_labels=no_labels,
        )
        self.labels = tree['ivl']
        self.figure.ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        self.figure.ax.get_xaxis().set_major_formatter(FuncFormatter(lambda x, p: format(int(x), ',')))

        xmax, xmin = self.figure.ax.set_xlim(xlim)

        if show_node_info:
            icoord, dcoord = tree['icoord'], tree['dcoord']
            for i, d in zip(icoord, dcoord):
                x = 0.5 * sum(i[1:3])
                y = d[1]
                if y <= xmax:
                    text = f"{y:.0f}"
                    self.figure.annotate(text, (y, x))
        if labels_color:
            ylbls = self.figure.ax.get_ymajorticklabels()
            for ylbl in ylbls:
                text = ylbl.get_text()
                color = labels_color.get(text, "#000000")
                ylbl.set_color(color)


def make_newick(node, newick, parentdist, leaf_names):
    """Convert scipy dendrogram to newick format."""
    if node.is_leaf():
        return "{}:{:.2f}{}".format(leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):{:.2f}{}".format(parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = make_newick(node.get_left(), newick, node.dist, leaf_names)
        newick = make_newick(node.get_right(), ",{}".format(newick), node.dist, leaf_names)
        newick = "({}".format(newick)
        return newick


def main():
    parser = argparse.ArgumentParser('cgMLST profiles cluster with single linkage algorithm.')
    parser.add_argument("-i", "--input", nargs='+',
                        required=True,
                        help="Path of cgMLST profile(s).")
    parser.add_argument("-o", "--output_path",
                        required=True,
                        help="Path of output directory.")
    args = parser.parse_args()

    png_filename = os.path.join(args.output_path, 'dendrogram.png')
    pdf_filename = os.path.join(args.output_path, 'dendrogram.pdf')
    svg_filename = os.path.join(args.output_path, 'dendrogram.svg')
    newick_filename = os.path.join(args.output_path, 'dendrogram.newick')
    dfs = []
    for file in args.input:
        filename = os.path.splitext(os.path.basename(file))[0]
        df = pd.read_csv(file, sep='\t', index_col=0, header=0, usecols=[0, 1], names=['locus_id', filename])
        dfs.append(df)
    profile = pd.concat(dfs, axis=1)
    dendrogram = Dendrogram(profile)
    dendrogram(show_node_info=True)
    for filename in (png_filename, pdf_filename, svg_filename):
        dendrogram.savefig(filename)
    dendrogram.to_newick(newick_filename)


if __name__ == '__main__':
    main()
