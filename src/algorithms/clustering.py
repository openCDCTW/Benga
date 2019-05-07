import fastcluster
import pandas as pd
import numpy as np
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


class DistanceMatrix:

    def __init__(self, profile):
        profile = profile.fillna('0')
        profile = profile.transpose()
        self._profile_values = profile.values
        self._profile_index = np.asarray(profile.index)
        self._distances = None

    @property
    def index(self):
        return self._profile_index

    @property
    def distance(self):
        if not self._distances:
            self._distances = pd.DataFrame()
            for index_1, item_1 in enumerate(self._profile_index, 0):
                for index_2, item_2 in enumerate(self._profile_index[index_1::], index_1):
                    self._distances.loc[item_1, item_2] = self._distances.loc[item_2, item_1] = hamming(
                        self._profile_values[index_1], self._profile_values[index_2]
                    )
        return self._distances


class Dendrogram:

    def __init__(self, dm, link):
        self._nodes = list(dm.index)
        self._newick = None
        if link == "single":
            self._linkage = fastcluster.single(squareform(dm.distance))
        elif link == "average":
            self._linkage = fastcluster.average(squareform(dm.distance))
        else:
            raise AttributeError("Invalid value {} for link in Dendrogram.".format(link))
        self._tree = hierarchy.to_tree(self._linkage, False)

    @property
    def newick(self):
        if not self._newick:
            self._newick = make_newick(self._tree, "", self._tree.dist, self._nodes)
        return self._newick

    def scipy_tree(self, file, distance_annotate=True, w=8, dpi=300):
        plt.style.use("fast")
        matplotlib.rcParams['lines.linewidth'] = 0.5
        fig, ax = plt.subplots(1, 1, figsize=(w, int(len(self._nodes)*0.3)))
        ax.grid(False)
        ax.patch.set_facecolor('none')
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        
        plt.rcParams['svg.fonttype'] = 'none'
        tree = hierarchy.dendrogram(self._linkage, labels=self._nodes, orientation="left",
                                    leaf_font_size=10, above_threshold_color="#000000", color_threshold=0)
        if distance_annotate:
            for i, d, in zip(tree['icoord'], tree['dcoord']):
                x = 0.5 * sum(i[1:3])
                y = d[1]
                plt.annotate(int(y), (y, x), xytext=(-2, 8), textcoords='offset points',
                             va='top', ha='right', fontsize=8)
        fig.savefig(file, dpi=dpi, bbox_inches='tight', pad_inches=1)

    def to_newick(self, file):
        with open(file, "w") as file:
            file.write(self.newick)


def hamming(item_1, item_2):
    return (item_1 != item_2).sum()


def make_newick(node, newick, parentdist, leaf_names):
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

