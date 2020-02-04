import fastcluster
import pandas as pd
import numpy as np
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from matplotlib import rcParams
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator, FuncFormatter

from ..utils.data import integer_encoding


class DistanceMatrix:
    def __init__(self, profile):
        assert isinstance(profile, pd.DataFrame)
        self.profile = profile
        shape = (len(profile.columns), len(profile.columns))
        self.matrix = np.empty(shape)

    def calculate(self):
        profile = self.profile.apply(integer_encoding, axis=1)
        values = profile.T.values
        for i_1, value_1 in enumerate(values):
            for i_2, value_2 in enumerate(values[i_1:], i_1):
                self.matrix[i_1, i_2] = self.matrix[i_2, i_1] = (value_1 != value_2).sum()
        return self.matrix


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


class Canvas:
    plt.style.use("fast")
    rcParams["lines.linewidth"] = 0.5

    def __init__(self, width, height):
        self.width = width
        self.height = height
        self.fig = None
        self.ax = None

    def _remove_edge(self):
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['bottom'].set_visible(False)
        self.ax.spines['left'].set_visible(False)
        self.ax.spines['right'].set_visible(False)

    def annotate(self, text, position, fontsize=8):
        self.ax.annotate(text, position, xytext=(-2, 8), textcoords='offset points', va='top', ha='right',
                         fontsize=fontsize)

    def draw(self):
        self.fig, self.ax = plt.subplots(1, 1, figsize=(self.width, self.height))
        self._remove_edge()
        self.ax.grid(False)
        self.ax.patch.set_facecolor('none')
        plt.close()

    def savefig(self, file, dpi=300):
        self.fig.savefig(file, dpi=dpi, bbox_inches='tight', pad_inches=1)


class Dendrogram:
    show_format = {'single': lambda x: '{:.0f}'.format(x),
                   'average': lambda x: '{:.1f}'.format(x)}

    def __init__(self, profile, linkage_method='single'):
        assert isinstance(profile, pd.DataFrame)
        self.profile = profile
        self.linkage_method = linkage_method
        self._linkage = None
        self.figure = None

    def to_newick(self, file):
        """Generate newick format with dendrogram."""
        if self._linkage is None:
            distmatrix = DistanceMatrix(self.profile).calculate()
            self._linkage = Linkage(distmatrix, self.linkage_method)
        tree = hierarchy.to_tree(self._linkage.matrix, False)
        newick = make_newick(tree, "", tree.dist, self.profile.columns)
        with open(file, 'w') as f:
            f.write(newick)

    def cluster(self, no_labels=False, show_node_info=False):
        """Generate dendrogram."""
        distmatrix = DistanceMatrix(self.profile).calculate()
        linkage = Linkage(distmatrix, self.linkage_method)
        canvas = Canvas(8, len(self.profile.columns, ) * 0.3)
        canvas.draw()
        dendrogram = hierarchy.dendrogram(
            linkage.matrix,
            ax=canvas.ax,
            labels=self.profile.columns,
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
                info = self.show_format[self.linkage_method](y)
                canvas.annotate(info, (y, x))
        canvas.ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        canvas.ax.get_xaxis().set_major_formatter(FuncFormatter(lambda x, p: format(int(x), ',')))
        self.figure = canvas

    def savefig(self, file):
        self.figure.savefig(file)


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

