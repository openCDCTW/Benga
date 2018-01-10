import fastcluster
import pandas as pd
import re
from ete3 import Tree
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')


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

    def make_tree(self, profile_file, names=None):
        profiles = pd.read_csv(profile_file, sep="\t", index_col=0)
        new_names = {}        
        if names:
            for i in names:
                new_names[i] = re.sub(r"_S\d+_L\d{3}[\d\w_-]+", "", names[i]).replace("-", ".")
        profiles.columns = list(map(lambda x: new_names[x], profiles.columns))
        self._nodes = list(profiles.columns)
        distances = distance_matrix(profiles)
        self._linkage = fastcluster.average(squareform(distances))
        self._tree = hierarchy.to_tree(self._linkage, False)

    def to_newick(self, file):
        with open(file, "w") as file:
            file.write(self.newick)


def hamming(xs, ys):
    return sum(xs.ne(ys) & ~(xs.isnull() & ys.isnull()))


def distance_matrix(profiles):
    distances = pd.DataFrame(index=profiles.columns, columns=profiles.columns)
    for x in profiles.columns:
        for y in profiles.columns:
            distances.loc[x, y] = hamming(profiles[x], profiles[y])
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

