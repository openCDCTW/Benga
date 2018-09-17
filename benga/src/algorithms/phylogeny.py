import fastcluster
import pandas as pd
from ete3 import Tree
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import Counter


class Dendrogram:

    def __init__(self, tree=None):
        self._nodes = None
        self._tree = tree
        self._newick = None
        self._ete_tree = None
        self._linkage = None

    @property
    def newick(self):
        if not self._newick:
            self._newick = make_newick(self._tree, "", self._tree.dist, self._nodes)
        return self._newick

    @property
    def ete_tree(self):
        if not self._ete_tree:
            self._ete_tree = Tree(self.newick)
        return self._ete_tree

    def render_on(self, file, w=900, h=1200, units="px", dpi=300, *args):
        self.ete_tree.render(file, w=w, h=h, units=units, dpi=dpi, *args)

    def scipy_tree(self, file, w=8, dpi=300):
        plt.style.use("ggplot")
        fig, ax = plt.subplots(1, 1, figsize=(w, int(len(self._nodes)*0.3)))
        ax.grid(False)
        ax.tick_params(axis='x', bottom='off', top='off', labelbottom='off')
        ax.patch.set_facecolor('none')
        plt.rcParams['svg.fonttype'] = 'none'
        hierarchy.dendrogram(self._linkage, labels=self._nodes, orientation="left",
                             leaf_font_size=10, above_threshold_color="#808080")
        fig.savefig(file, dpi=dpi, bbox_inches='tight', pad_inches=1)

    def make_tree(self, profile_file):
        profiles = pd.read_csv(profile_file, sep="\t", index_col=0)
        self._nodes = list(profiles.columns)
        distances = distance_matrix(profiles)
        self._linkage = fastcluster.average(squareform(distances))
        self._tree = hierarchy.to_tree(self._linkage, False)

    def to_newick(self, file):
        with open(file, "w") as file:
            file.write(self.newick)


def hamming(xs, ys):
    return sum(xs.ne(ys) & ~(xs.isnull() & ys.isnull()))


def distance(array_1, array_2, name_1, name_2):
    d = Counter(array_1 == array_2)[False]
    return name_1, name_2, d


def distance_matrix_new(profile):
    profile = profile.fillna('0')
    profile = profile.transpose()
    profile_array = profile.values
    profile_index = [i for i in range(len(profile.index))]
    pairs = []
    for i in range(len(profile_index)):
        index_1 = profile_index.pop(0)
        pairs.append((profile.index[index_1], profile.index[index_1], 0))
        for index_2 in profile_index:
            dis = distance(profile_array[index_1], profile_array[index_2], profile.index[index_1], profile.index[index_2])
            pairs.append(dis)
    distances = pd.DataFrame()
    for pair in pairs:
        distances.loc[pair[0], pair[1]] = distances.loc[pair[1], pair[0]] = pair[2]
    return distances


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

