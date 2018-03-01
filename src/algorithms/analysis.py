import json
import os.path
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from src.utils import db
plt.style.use("ggplot")


TO_INT = ["No. isolates", "No. sequences", "Genome Fragment", "Order within Fragment", "Accessory Fragment",
          "Accessory Order with Fragment", "Min group size nuc", "Max group size nuc", "Avg group size nuc"]
TO_STR = ["Non-unique Gene name", "Annotation", "QC"]
TO_FLOAT = ["Avg sequences per isolate"]


def power():
    counts = db.from_sql("select locus_id, count(locus_id) as counts from pairs group by locus_id;")
    counts["log_counts"] = np.log2(counts["counts"])
    return np.sum(counts["log_counts"])


def locus_entropy(x):
    prob = x / np.sum(x)
    return np.sum(prob * np.log2(prob))


def richness(weighted=True):
    sql = "select a.locus_id, a.allele_id, b.count" \
          " from pairs as a" \
          " left join (select allele_id, count from alleles) as b" \
          " on a.allele_id=b.allele_id;"
    counts = db.from_sql(sql)
    ent = counts.groupby("locus_id").agg({"count": locus_entropy})
    if weighted:
        sql = "select locus_id, occurrence from loci;"
        loci = db.from_sql(sql)
        weight = pd.merge(ent, loci, left_index=True, right_on="locus_id")
        return np.average(weight["count"], weights=weight["occurrence"])
    else:
        return np.average(ent)


def calculate_loci_coverage(input_dir, output_dir, database):
    db.load_database_config()
    subject_number = count_subjects(input_dir)
    plot_stats(output_dir, subject_number, database)


def count_subjects(input_dir):
    input_file = os.path.join(input_dir, "allele_profiles.tsv")
    with open(input_file, "r") as file:
        first_line = file.readline()
    genomes = first_line.strip().split("\t")
    return len(genomes)


def plot_stats(output_dir, subject_number, database):
    sql = "select locus_id, num_isolates, is_paralog from locus_meta where is_paralog=TRUE;"
    table = db.from_sql(sql, database=database)
    table["owned by"] = [int(x / subject_number * 100) for x in table["num_isolates"]]
    plot_genome_coverage(table["owned by"], output_dir, perc=0)
    plot_genome_coverage(table["owned by"], output_dir)
    plot_genome_coverage(table["owned by"], output_dir, perc=0, cumulative=-1)
    plot_genome_coverage(table["owned by"], output_dir, cumulative=-1)


def plot_genome_coverage(data, output_dir, perc=5, cumulative=False):
    prefix = "cumulative_genome_coverage" if cumulative else "genome_coverage"
    pic_name = "{}_{}_prec.png".format(prefix, perc) if perc != 0 else "{}.png".format(prefix)
    output_file = os.path.join(output_dir, pic_name)
    title = "Cumulative genome coverage distribution" if cumulative else "Genome coverage distribution"

    # plot
    fig = plt.figure(figsize=(12, 9))
    plt.hist(data[data >= perc], bins=50, cumulative=cumulative, histtype="step", lw=2)
    plt.title(title, fontsize=25)
    plt.xlabel("Percentage of genomes covered by loci (%)", fontsize=18)
    plt.ylabel("Number of locus", fontsize=18)
    fig.savefig(output_file)


# allele


def calculate_allele_length(input_dir, output_dir):
    db.load_database_config()
    db_dir = os.path.join(input_dir, "database")
    plot_length_heamap(db_dir, output_dir)


def plot_length_heamap(input_dir, output_dir, interval=20):
    output_file = os.path.join(output_dir, "allele_length_heatmap.png")
    allele_freq = parse_freq(input_dir)
    allele_len = parse_length(os.path.join(input_dir, "locusfiles"), allele_freq)
    table = to_dataframe(allele_len, interval)
    table = table.apply(lambda x: 100 * x / np.sum(x), axis=1)

    # sort by scheme order
    freq = pd.read_csv(os.path.join(input_dir, "scheme.tsv"), sep="\t", usecols=["locus"])
    table = pd.merge(freq, table, left_on="locus", right_index=True).set_index("locus")

    table = table.apply(mask_by_length, axis=1).apply(np.floor, axis=1)
    to_show = table.iloc[0:100, 0:80]

    # plot
    fig = plt.figure(figsize=(24, 16))
    ax = sns.heatmap(to_show, annot=True, annot_kws={}, mask=to_show.isnull(),
                     cmap=sns.blend_palette(["#446e8c", "#f6ff6d"], n_colors=20, as_cmap=True))
    fig.add_axes(ax)
    plt.xticks(rotation=30)
    plt.yticks(rotation=0)
    plt.savefig(output_file)


def to_dataframe(allele_len, interval):
    collect = []
    for locus, lens in allele_len.items():
        times = [(int(x / interval) + 1) * interval for x in lens]
        series = pd.Series(data=Counter(times), name=locus)
        collect.append(series)
    table = pd.concat(collect, axis=1).fillna(0).T
    return table


def parse_freq(db_dir):
    file = os.path.join(db_dir, "allele_frequency.json")
    with open(file, "r") as f:
        allele_freq = json.loads(f.read())
    return allele_freq


def parse_length(locus_dir, freq):
    allele_len = {}
    for fastafile in os.listdir(locus_dir):
        file = os.path.join(locus_dir, fastafile)
        locus = fastafile.split(".")[0]
        allele_f = freq[locus]
        weighted = []
        for rec in SeqIO.parse(file, "fasta"):
            l = len(rec.seq)
            weighted.extend([l] * allele_f[rec.id])
        allele_len[locus] = weighted
    return allele_len


def mask_by_length(x):
    y = x.copy()
    tostop = 0
    for i in reversed(y.index):
        if y[i] != 0.0:
            tostop = i
            break
    y[tostop+1:] = np.nan
    return y
