import os.path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from src.models import logs
from src.utils import db
plt.style.use("ggplot")


TO_INT = ["No. isolates", "No. sequences", "Genome Fragment", "Order within Fragment", "Accessory Fragment",
          "Accessory Order with Fragment", "Min group size nuc", "Max group size nuc", "Avg group size nuc"]
TO_STR = ["Non-unique Gene name", "Annotation", "QC"]
TO_FLOAT = ["Avg sequences per isolate"]


def power(database):
    lf = logs.LoggerFactory()
    lf.addConsoleHandler()
    logger = lf.create()
    db.load_database_config(logger=logger)
    sql = "select locus_id, count(locus_id) as counts from pairs group by locus_id;"
    counts = db.from_sql(sql, database=database)
    counts["log_counts"] = np.log2(counts["counts"])
    return np.sum(counts["log_counts"])


def locus_entropy(x):
    prob = x / np.sum(x)
    return np.sum(prob * np.log2(prob))


def richness(database, weighted=True):
    lf = logs.LoggerFactory()
    lf.addConsoleHandler()
    logger = lf.create()
    db.load_database_config(logger=logger)
    sql = "select a.locus_id, a.allele_id, b.count" \
          " from pairs as a" \
          " left join (select allele_id, count from alleles) as b" \
          " on a.allele_id=b.allele_id;"
    counts = db.from_sql(sql, database=database)
    ent = counts.groupby("locus_id").agg({"count": locus_entropy})
    if weighted:
        sql = "select locus_id, occurrence from loci;"
        loci = db.from_sql(sql, database=database)
        weight = pd.merge(ent, loci, left_index=True, right_on="locus_id")
        return np.average(weight["count"], weights=weight["occurrence"])
    else:
        return np.average(ent)


def calculate_loci_coverage(input_dir, output_dir, database):
    lf = logs.LoggerFactory()
    lf.addConsoleHandler()
    logger = lf.create()
    db.load_database_config(logger=logger)
    subject_number = count_subjects(input_dir)
    plot_stats(output_dir, subject_number, database)


def count_subjects(input_dir):
    input_file = os.path.join(input_dir, "allele_profiles.tsv")
    with open(input_file, "r") as file:
        first_line = file.readline()
    genomes = first_line.strip().split("\t")
    return len(genomes)


def plot_stats(output_dir, subject_number, database):
    sql = "select locus_id, num_isolates, is_paralog from locus_meta where is_paralog=FALSE;"
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


def calculate_allele_length(output_dir, database, interval=20):
    lf = logs.LoggerFactory()
    lf.addConsoleHandler()
    logger = lf.create()
    db.load_database_config(logger=logger)
    plot_length_heamap(output_dir, database, interval=interval)


def plot_length_heamap(output_dir, database, interval):
    output_file = os.path.join(output_dir, "allele_length_heatmap.png")
    allele_info = get_allele_info(database)
    allele_info["intervals"] = list(map(lambda x: (int(x / interval) + 1) * interval, allele_info["length"]))
    pairs = db.from_sql("select * from pairs;", database=database)
    collect = []
    for locus_id, df in pairs.groupby("locus_id"):
        df2 = pd.merge(df, allele_info, on="allele_id", how="left")
        series = df2.groupby("intervals")["count"].sum()
        series.name = locus_id
        collect.append(series)
    table = pd.concat(collect, axis=1).fillna(0).T
    table = table.apply(lambda x: 100 * x / np.sum(x), axis=1)

    # sort by scheme order
    freq = db.from_sql("select locus_id from loci order by occurrence DESC;", database=database)
    table = pd.merge(freq, table, left_on="locus_id", right_index=True).set_index("locus_id")

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


def get_allele_info(database):
    sql = "select allele_id, char_length(dna_seq) as length, count from alleles;"
    return db.from_sql(sql, database=database)


def mask_by_length(x):
    y = x.copy()
    tostop = 0
    for i in reversed(y.index):
        if y[i] != 0.0:
            tostop = i
            break
    y[tostop+1:] = np.nan
    return y
