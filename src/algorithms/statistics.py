import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from src.models import logs
from src.utils import db, seq, files
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
    logger.info("Start calculating locus coverage...")
    subject_number = count_subjects(input_dir)
    logger.info("Start plotting locus coverage...")
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
    logger.info("Start calculating allele length heatmap...")
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


def reference_self_blastp(output_dir, database):
    query = "select loci.locus_id, alleles.peptide_seq" \
            " from loci inner join alleles on loci.ref_allele=alleles.allele_id;"
    ref = db.from_sql(query, database=database)
    ref_recs = [seq.new_record(row["locus_id"], row["peptide_seq"]) for _, row in ref.iterrows()]
    ref_faa = files.joinpath(output_dir, "ref_seq.faa")
    seq.save_records(ref_recs, ref_faa)

    ref_db = files.joinpath(output_dir, "ref_db")
    seq.compile_blastpdb(ref_faa, ref_db)

    blastp_out_file = files.joinpath(output_dir, "ref_db.blastp.out")
    seq.query_blastpdb(ref_faa, ref_db, blastp_out_file, seq.BLAST_COLUMNS)
    return blastp_out_file


def identify_pairs(df):
    sseqids = df["sseqid"].tolist()
    pairs = []
    counted = set()
    for _, row in df.iterrows():
        if row["qseqid"] in sseqids and not row["qseqid"] in counted:
            counted.update(row["qseqid"])
            counted.update(row["sseqid"])
            pairs.append((row["qseqid"], row["sseqid"]))
    return pairs


def collect_high_occurrence_loci(pairs, database):
    query = "select locus_id, occurrence from loci;"
    occur = db.from_sql(query, database=database)
    drops = set()
    for id1, id2 in pairs:
        ocr1 = occur.loc[occur["locus_id"] == id1, "occurrence"].iloc[0]
        ocr2 = occur.loc[occur["locus_id"] == id2, "occurrence"].iloc[0]
        drops.update([id2] if ocr1 >= ocr2 else [id1])
    filtered_loci = set(occur["locus_id"]) - drops
    return filtered_loci


def query_covs(covs, key):
    if key in covs.keys():
        return covs[key]
    else:
        return np.nan


def filter_locus(blastp_out_file, database):
    blastp_out = pd.read_csv(blastp_out_file, sep="\t", header=None, names=seq.BLAST_COLUMNS)
    blastp_out = blastp_out[blastp_out["pident"] >= 95]
    blastp_out = blastp_out[blastp_out["qseqid"] != blastp_out["sseqid"]]
    covs = {(row["qseqid"], row["sseqid"]): row["qcovs"] for _, row in blastp_out.iterrows()}
    blastp_out["scovs"] = [query_covs(covs, (row["sseqid"], row["qseqid"])) for _, row in blastp_out.iterrows()]
    blastp_out["qlen/alen"] = 100 / blastp_out["qcovs"]
    blastp_out["qlen/slen"] = blastp_out["scovs"] / blastp_out["qcovs"]
    blastp_out = blastp_out[(0.75 < blastp_out["qlen/alen"]) & (blastp_out["qlen/alen"] <= 1.25)]
    blastp_out = blastp_out[(0.75 < blastp_out["qlen/slen"]) & (blastp_out["qlen/slen"] <= 1.25)]
    pairs = identify_pairs(blastp_out)
    filtered_loci = collect_high_occurrence_loci(pairs, database)
    return filtered_loci


def extract_database_ref_sequence(old_database, new_database, keep_loci):
    alleles = db.from_sql("select * from alleles;", database=old_database)
    loci = db.from_sql("select * from loci;", database=old_database)
    locus_meta = db.from_sql("select * from locus_meta;", database=old_database)
    pairs = db.from_sql("select * from pairs;", database=old_database)

    ref_alleles = set(loci["ref_allele"])
    alleles = alleles[alleles["allele_id"].isin(ref_alleles)]
    pairs = pairs[(pairs["allele_id"].isin(ref_alleles)) & (pairs["locus_id"].isin(keep_loci))]
    loci = loci[loci["locus_id"].isin(keep_loci)]

    db.append_to_sql("alleles", alleles, new_database)
    db.append_to_sql("locus_meta", locus_meta, new_database)
    db.append_to_sql("loci", loci, new_database)
    db.append_to_sql("pairs", pairs, new_database)


def build_locus_library(output_dir, old_database, logger=None):
    if not logger:
        lf = logs.LoggerFactory()
        lf.addConsoleHandler()
        lf.addFileHandler(files.joinpath(output_dir, "make_database.log"))
        logger = lf.create()
    db.load_database_config(logger=logger)

    new_database = old_database + "_new"
    db.createdb(new_database)
    db.create_pgadb_relations(new_database)

    blastp_out_file = reference_self_blastp(output_dir, old_database)
    filtered_loci = filter_locus(blastp_out_file, old_database)
    extract_database_ref_sequence(old_database, new_database, filtered_loci)
    os.remove(blastp_out_file)
    return new_database
