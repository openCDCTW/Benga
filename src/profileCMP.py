import shutil
import os
import pandas as pd
import matplotlib.pyplot as plt
import fastcluster
from scipy.cluster.hierarchy import dendrogram
from src.utils import *

dest_dir = "/media/pika/Workbench/workspace/forYH_MSGAdb-builder/04--ProfileCMP"
src_dir = "/media/pika/Workbench/workspace/forYH_MSGAdb-builder/03--wgProfileing/test"
scheme_selection = "locusAP_core_100p_dispensable"

# cpScheme

shutil.copy(joinpath(src_dir, "files.new.list"), ".")
shutil.copytree(joinpath(src_dir, "OUTPUT"), "scheme")



# makeDiffMatrix

scheme_path = {"locusAP_core": "/scheme/locusAP/locusAP_core",
               "locusAP_core_10p_dispensable": "/scheme/locusAP/locusAP_core_10p_dispensable",
               "locusAP_core_20p_dispensable": "/scheme/locusAP/locusAP_core_20p_dispensable",
               "locusAP_core_30p_dispensable": "/scheme/locusAP/locusAP_core_30p_dispensable",
               "locusAP_core_40p_dispensable": "/scheme/locusAP/locusAP_core_40p_dispensable",
               "locusAP_core_50p_dispensable": "/scheme/locusAP/locusAP_core_50p_dispensable",
               "locusAP_core_60p_dispensable": "/scheme/locusAP/locusAP_core_60p_dispensable",
               "locusAP_core_70p_dispensable": "/scheme/locusAP/locusAP_core_70p_dispensable",
               "locusAP_core_80p_dispensable": "/scheme/locusAP/locusAP_core_80p_dispensable",
               "locusAP_core_90p_dispensable": "/scheme/locusAP/locusAP_core_90p_dispensable",
               "locusAP_core_100p_dispensable": "/scheme/locusAP/locusAP_core_100p_dispensable",
               "locusAP_pan": "/scheme/locusAP/locusAP_pan",
               "locusAP_unique": "/scheme/locusAP/locusAP_unique",
               "wgMLST_core": "/scheme/wgMLST/wgMLST_core",
               "wgMLST_core_10p_dispensable": "/scheme/wgMLST/wgMLST_core_10p_dispensable",
               "wgMLST_core_20p_dispensable": "/scheme/wgMLST/wgMLST_core_20p_dispensable",
               "wgMLST_core_30p_dispensable": "/scheme/wgMLST/wgMLST_core_30p_dispensable",
               "wgMLST_core_40p_dispensable": "/scheme/wgMLST/wgMLST_core_40p_dispensable",
               "wgMLST_core_50p_dispensable": "/scheme/wgMLST/wgMLST_core_50p_dispensable",
               "wgMLST_core_60p_dispensable": "/scheme/wgMLST/wgMLST_core_60p_dispensable",
               "wgMLST_core_70p_dispensable": "/scheme/wgMLST/wgMLST_core_70p_dispensable",
               "wgMLST_core_80p_dispensable": "/scheme/wgMLST/wgMLST_core_80p_dispensable",
               "wgMLST_core_90p_dispensable": "/scheme/wgMLST/wgMLST_core_90p_dispensable",
               "wgMLST_core_100p_dispensable": "/scheme/wgMLST/wgMLST_core_100p_dispensable",
               "wgMLST_pan": "/scheme/wgMLST/wgMLST_pan",
               "wgMLST_unique": "/scheme/wgMLST/wgMLST_unique"}

target_dir = dest_dir + scheme_path[scheme_selection]

namemap = {}
for file in open("files.new.list", "r").read().splitlines():
    a1, a2 = file.split("\t")[0:2]
    newname = a1.split(".")[0]
    oldname = a2.split(".")[0]
    namemap[newname] = oldname


# collect labels
nodeLabel = []
for filename in open(joinpath(target_dir, "files.list"), "r").read().splitlines():
    fn = filename.split(".")[1]
    nodeLabel.append(fn.replace("Assembly_", "A") + "\t" + namemap[fn])

with open("nodeLableMatch.txt", "w") as file:
    file.write("\n".join(nodeLabel))


# collect profiles
matrix = pd.DataFrame()
for filename in open(joinpath(target_dir, "files.list"), "r").read().splitlines():
    column = pd.Series()
    for line in open(target_dir + "/" + filename).read().splitlines():
        if line == "":
            continue
        locus, hit = line.split("\t")
        column[locus] = True if hit != 0 else False
    sample = filename.split(".")[1]
    matrix[sample] = column

matrix.to_csv("profiles.tsv", sep="\t")


Y = fastcluster.linkage(matrix, method="average", metric="hamming")
dendrogram(Y)
plt.show()
plt.close()

