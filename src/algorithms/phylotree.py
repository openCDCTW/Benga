import fastcluster
import pandas as pd
from ete3 import Tree
from scipy.cluster.hierarchy import to_tree

from src.utils import files


def linkage2newick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "{}:{:.2f}{}".format(leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):{:.2f}{}".format(parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = linkage2newick(node.get_left(), newick, node.dist, leaf_names)
        newick = linkage2newick(node.get_right(), ",{}".format(newick), node.dist, leaf_names)
        newick = "({}".format(newick)
        return newick


def hamming(xs, ys):
    return sum(x != y for x, y in zip(xs, ys))


def make_tree(input_dir, scheme_selection):
    schemes = ["wgMLST_core", "wgMLST_pan"]
    scheme_path = {s: files.joinpath("scheme", s.split("_")[0], s) for s in schemes}
    target_file = files.joinpath(input_dir, scheme_path[scheme_selection] + ".tsv")
    profiles = pd.read_csv(target_file, sep="\t", index_col=0)

    ## build phylogenetic tree
    # distance matrix using hamming distance
    distances = pd.DataFrame(index=profiles.columns, columns=profiles.columns)
    for x in profiles.columns:
        for y in profiles.columns:
            distances.loc[x, y] = hamming(profiles[x], profiles[y])

    # linkage
    linkage = fastcluster.average(distances)
    tree = to_tree(linkage, False)
    newick = linkage2newick(tree, "", tree.dist, list(profiles.columns))
    return Tree(newick), newick


def to_graph(tree, file):
    tree.render(file, w=800, h=600, units="px", dpi=300)


def to_newick(newick, file):
    with open(file, "w") as file:
        file.write(newick)
