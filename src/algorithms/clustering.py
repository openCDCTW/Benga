import fastcluster
import pandas as pd
import numpy as np
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from matplotlib import rcParams
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator, FuncFormatter

from ..utils.data import integer_encoding


class Distance:
    def __init__(self, profile):
        assert isinstance(profile, pd.DataFrame)
        self.profile = profile.apply(integer_encoding, axis=1)
        shape = (len(profile.columns), len(profile.columns))
        self.matrix = np.empty(shape)

    def calculate(self):
        values = self.profile.T.values
        for i_1, value_1 in enumerate(values):
            for i_2, value_2 in enumerate(values[i_1:], i_1):
                self.matrix[i_1, i_2] = self.matrix[i_2, i_1] = (value_1 != value_2).sum()


class Linkage:
    def __init__(self, distmatrix, method='single'):
        self.cdm = squareform(distmatrix)
        self.method = method
        if method == 'single':
            self.matrix = fastcluster.single(self.cdm)
        elif method == 'average':
            self.matrix = fastcluster.average(self.cdm)
        else:
            raise


class Figure:
    plt.style.use("fast")
    rcParams["lines.linewidth"] = 0.5

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

    def savefig(self, file, dpi=300):
        self.fig.savefig(file, dpi=dpi, bbox_inches='tight', pad_inches=1)


class Dendrogram:
    show_format = {'single': lambda x: '{:.0f}'.format(x),
                   'average': lambda x: '{:.1f}'.format(x)}

    def __init__(self, profile, linkage_method='single'):
        assert isinstance(profile, pd.DataFrame)
        self.distance = Distance(profile)
        self.distance.calculate()
        self.linkage = Linkage(distmatrix=self.distance.matrix, method=linkage_method)
        self.figure = Figure(12, len(self.distance.profile.columns) * 0.3)

    def to_newick(self, file):
        """Generate newick format with dendrogram."""
        tree = hierarchy.to_tree(self.linkage.matrix, False)
        newick = make_newick(tree, "", tree.dist, self.distance.profile.columns)
        with open(file, 'w') as f:
            f.write(newick)

    def cluster(self, no_labels=False, show_node_info=False):
        """Generate dendrogram."""
        dendrogram = hierarchy.dendrogram(
            self.linkage.matrix,
            ax=self.figure.ax,
            labels=self.distance.profile.columns,
            orientation="left",
            leaf_font_size=12,
            above_threshold_color="#000000",
            color_threshold=0,
            no_labels=no_labels,
        )
        if show_node_info:
            icoord, dcoord = dendrogram['icoord'], dendrogram['dcoord']
            for i, d in zip(icoord, dcoord):
                x = 0.5 * sum(i[1:3])
                y = d[1]
                info = self.show_format[self.linkage.method](y)
                self.figure.annotate(info, (y, x))
        self.figure.ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        self.figure.ax.get_xaxis().set_major_formatter(FuncFormatter(lambda x, p: format(int(x), ',')))


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

