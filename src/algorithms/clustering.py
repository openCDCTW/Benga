import fastcluster
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import Counter
from matplotlib.ticker import MaxNLocator


class DistanceMatrix:

    def __init__(self, profile):
        profile = profile.fillna('0')
        profile = profile.transpose()
        self._profile_array = profile.values
        self._profile_index = list(profile.index)
        self._distances = None

    @property
    def index(self):
        return self._profile_index

    @property
    def distance(self):
        if not self._distances:
            pairs = []
            for i in range(len(self._profile_index)):
                data = (self._profile_index[i], self._profile_index[i], 0)
                pairs.append(data)
                for j in range(i+1, len(self._profile_index)):
                    dist = hamming(self._profile_array[i], self._profile_array[j])
                    data = (self._profile_index[i], self._profile_index[j], dist)
                    pairs.append(data)
            self._distances = pd.DataFrame()
            for pair in pairs:
                self._distances.loc[pair[0], pair[1]] = self._distances.loc[pair[1], pair[0]] = pair[2]
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


def hamming(array_1, array_2):
    return Counter(array_1 == array_2)[False]


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

