import json
import os
import re
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor

import pandas as pd
from Bio import SeqIO
from benga.src.utils import seq, files, cmds, operations, db, logs


def move_file(annotate_dir, dest_dir, ext):
    for folder in os.listdir(annotate_dir):
        path = files.joinpath(annotate_dir, folder)
        for file in os.listdir(path):
            current_file = files.joinpath(path, file)
            if file.endswith(ext):
                os.rename(current_file, files.joinpath(dest_dir, file))


def create_noncds(database_dir, gff_dir):
    noncds = defaultdict(list)
    for file in os.listdir(gff_dir):
        name, ext = file.split(".")
        for line in open(files.joinpath(gff_dir, file), "r").read().splitlines():
            if line.startswith("##sequence") or line.startswith("##gff"):
                continue
            if line.startswith("##FASTA"):
                break

            token = line.split("\t")
            seq_type, annotation = token[2], token[8]

            if annotation.startswith("ID="):
                prokkaid = annotation.split(";")[0][3:]
                if seq_type != "CDS":
                    noncds[name].append(prokkaid)
    with open(files.joinpath(database_dir, "nonCDS.json"), "w") as file:
        file.write(json.dumps(noncds))


def filter_tRNA(matrix):
    c = re.compile(r"tRNA-\w+\(\w{3}\)")
    is_trna = []
    for x in matrix["description"]:
        m = c.match(x)
        is_trna.append(m is not None)
    matrix["is_trna"] = is_trna
    return matrix[~matrix["is_trna"]].drop("is_trna", axis=1)


def filter_rRNA(matrix):
    c1 = re.compile("ribosomal RNA")
    c2 = re.compile("subunit")
    is_rrna = []
    for x in matrix["description"]:
        m1 = c1.match(x)
        if m1 is not None:
            m2 = c2.match(x)
            is_rrna.append(m2 is None)
        else:
            is_rrna.append(False)
    matrix["is_rrna"] = is_rrna
    return matrix[~matrix["is_rrna"]].drop("is_rrna", axis=1)


def extract_profiles(roary_matrix_file, dbname, metadata_cols=13):
    matrix = pd.read_csv(roary_matrix_file)
    matrix["Gene"] = matrix["Gene"].str.replace("/", "_")
    matrix["Gene"] = matrix["Gene"].str.replace(" ", "_")
    rename_cols = {"Gene": "locus_id", "No. isolates": "num_isolates", "No. sequences": "num_sequences",
                   "Annotation": "description"}
    matrix.rename(columns=rename_cols, inplace=True)
    matrix.set_index("locus_id", inplace=True)
    matrix = filter_tRNA(matrix)
    matrix = filter_rRNA(matrix)
    save_locus_metadata(matrix, dbname)
    profiles = matrix.iloc[:, metadata_cols:]
    isolates = len(matrix.columns) - metadata_cols
    return profiles, isolates


def save_locus_metadata(matrix, dbname, select_col=None, repeat_tol=1.2):
    if not select_col:
        select_col = ["locus_id", "num_isolates", "num_sequences", "description", "is_paralog"]
    avg = "Avg sequences per isolate"
    meta = matrix.copy()
    meta["is_paralog"] = [x > repeat_tol for x in meta[avg]]
    meta = meta.reset_index()[select_col]
    db.append_to_sql("locus_meta", meta, dbname)


def collect_allele_infos(profiles, ffn_dir):
    freq = defaultdict(Counter)
    new_profiles = []
    for subject, profile in profiles.iteritems():
        ffn_file = files.joinpath(ffn_dir, "{}.ffn".format(subject))
        seqs = {record.id: record.seq for record in SeqIO.parse(ffn_file, "fasta")}

        new_profile = pd.Series(name=subject)
        for locus, prokka_str in profile.dropna().iteritems():
            if "\t" not in prokka_str:
                allele = seqs[prokka_str]
                freq[locus].update([allele])
                new_profile.at[locus] = operations.make_seqid(allele)
            else:
                prokka_id = prokka_str.split("\t")
                alleles = [seqs[x] for x in prokka_id]
                freq[locus].update(alleles)
                v = "\t".join(operations.make_seqid(x) for x in alleles)
                new_profile.at[locus] = v
        new_profiles.append(new_profile)
    new_profiles = pd.concat(new_profiles, axis=1).sort_index().sort_index(axis=1)
    return new_profiles, freq


def reference_self_blastp(output_dir, freq):
    ref_recs = [seq.new_record(locus, counter.most_common(1)[0][0].translate(table=11)) for locus, counter in freq.items()]
    ref_length = {rec.id: len(rec.seq) for rec in ref_recs}
    ref_faa = files.joinpath(output_dir, "ref_seq.faa")
    seq.save_records(ref_recs, ref_faa)

    ref_db = files.joinpath(output_dir, "ref_db")
    seq.compile_blastpdb(ref_faa, ref_db)

    blastp_out_file = files.joinpath(output_dir, "ref_db.blastp.out")
    seq.query_blastpdb(ref_faa, ref_db, blastp_out_file, seq.BLAST_COLUMNS)
    return blastp_out_file, ref_length

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


def select_drop_loci(df):
    drop1 = df[df["locus_id"].str.contains("group_")]["locus_id"]
    drop2 = df[df["locus_id"].str.contains(r"_\d$")]["locus_id"]
    return set(drop1) | set(drop2)


def collect_high_occurrence_loci(pairs, total_isolates, drop_by_occur):
    occur = db.from_sql("select locus_id, num_isolates from locus_meta;")
    occur["occurrence"] = list(map(lambda x: round(x / total_isolates * 100, 2), occur["num_isolates"]))
    drops = set()
    for id1, id2 in pairs:
        ocr1 = occur.loc[occur["locus_id"] == id1, "occurrence"].iloc[0]
        ocr2 = occur.loc[occur["locus_id"] == id2, "occurrence"].iloc[0]
        drops.update([id2] if ocr1 >= ocr2 else [id1])
    drops2 = select_drop_loci(occur[occur["occurrence"] < drop_by_occur])
    filtered_loci = set(occur["locus_id"]) - drops - drops2
    return filtered_loci


def filter_locus(blastp_out_file, ref_length, total_isolates, drop_by_occur):
    blastp_out = pd.read_csv(blastp_out_file, sep="\t", header=None, names=seq.BLAST_COLUMNS)
    blastp_out = blastp_out[blastp_out["pident"] >= 95]
    blastp_out = blastp_out[blastp_out["qseqid"] != blastp_out["sseqid"]]
    blastp_out["qlen"] = list(map(lambda x: ref_length[x], blastp_out["qseqid"]))
    blastp_out["slen"] = list(map(lambda x: ref_length[x], blastp_out["sseqid"]))
    blastp_out["qlen/slen"] = blastp_out["qlen"] / blastp_out["slen"]
    blastp_out["qlen/alen"] = blastp_out["qlen"] / blastp_out["length"]
    blastp_out = blastp_out[(0.75 <= blastp_out["qlen/slen"]) & (blastp_out["qlen/slen"] < 1.25)]
    blastp_out = blastp_out[(0.75 <= blastp_out["qlen/alen"]) & (blastp_out["qlen/alen"] < 1.25)]
    pairs = identify_pairs(blastp_out)
    filtered_loci = collect_high_occurrence_loci(pairs, total_isolates, drop_by_occur)
    return filtered_loci


def to_allele_table(data, dbname):
    df = pd.DataFrame(data, columns=["allele_id", "dna_seq", "peptide_seq", "count"])
    df = df.groupby("allele_id").agg({"dna_seq": "first", "peptide_seq": "first", "count": "sum"})
    df.reset_index(inplace=True)
    db.append_to_sql("alleles", df, dbname)


def to_pair_table(data, dbname):
    df = pd.DataFrame(data, columns=["allele_id", "locus_id"])
    db.append_to_sql("pairs", df, dbname)


def save_sequences(freq, dbname):
    alleles = []
    pairs = []
    for locus, counter in freq.items():
        for allele, count in counter.items():
            dna_seq = str(allele)
            pept_seq = str(allele.translate(table=11))
            allele_id = operations.make_seqid(dna_seq)
            alleles.append((allele_id, dna_seq, pept_seq, count))
            pairs.append((allele_id, locus))
    to_allele_table(alleles, dbname)
    to_pair_table(pairs, dbname)


def make_schemes(freq, total_isolates):
    refseqs = {locus: operations.make_seqid(counter.most_common(1)[0][0]) for locus, counter in freq.items()}
    schemes = db.from_sql("select locus_id, num_isolates from locus_meta;")
    schemes = schemes[schemes["locus_id"].isin(refseqs.keys())]
    schemes["occurrence"] = list(map(lambda x: round(x/total_isolates * 100, 2), schemes["num_isolates"]))
    schemes["ref_allele"] = list(map(lambda x: refseqs[x], schemes["locus_id"]))
    schemes = schemes[["locus_id", "occurrence", "ref_allele"]]
    db.append_to_sql("loci", schemes)


def annotate_configs(input_dir, output_dir, logger=None, threads=8):
    if not logger:
        lf = logs.LoggerFactory()
        lf.addConsoleHandler()
        lf.addFileHandler(files.joinpath(output_dir, "annotation.log"))
        logger = lf.create()

    logger.info("Formating contigs...")
    genome_dir = files.joinpath(output_dir, "Genomes")
    files.create_if_not_exist(genome_dir)
    contighandler = files.ContigHandler()
    contighandler.new_format(input_dir, genome_dir, replace_ext=False)
    namemap = contighandler.namemap
    with open(files.joinpath(output_dir, "namemap.json"), "w") as f:
        f.write(json.dumps(namemap))

    logger.info("Annotating...")
    annotate_dir = files.joinpath(output_dir, "Annotated")
    files.create_if_not_exist(annotate_dir)
    c = [cmds.form_prokka_cmd(x, genome_dir, annotate_dir) for x in namemap.keys()]
    with ProcessPoolExecutor(int(threads / 2)) as executor:
        executor.map(cmds.execute_cmd, c)

    logger.info("Moving protein CDS (.ffn) files...")
    ffn_dir = files.joinpath(output_dir, "FFN")
    files.create_if_not_exist(ffn_dir)
    move_file(annotate_dir, ffn_dir, ".ffn")

    logger.info("Moving annotation (.gff) files...")
    gff_dir = files.joinpath(output_dir, "GFF")
    files.create_if_not_exist(gff_dir)
    move_file(annotate_dir, gff_dir, ".gff")

    logger.info("Creating nonCDS.json...")
    create_noncds(output_dir, gff_dir)


def make_database(output_dir, drop_by_occur, logger=None, threads=2):
    if not logger:
        lf = logs.LoggerFactory()
        lf.addConsoleHandler()
        lf.addFileHandler(files.joinpath(output_dir, "make_database.log"))
        logger = lf.create()
    db.load_database_config(logger=logger)

    logger.info("Calculating the pan genome...")
    min_identity = 95
    c = cmds.form_roary_cmd(files.joinpath(output_dir, "GFF"), output_dir, min_identity, threads)
    logger.info("Run roary with following command: " + c)
    subprocess.run(c, shell=True)

    logger.info("Creating database")
    dbname = os.path.basename(output_dir[:-1] if output_dir.endswith("/") else output_dir)
    db.createdb(dbname)
    db.create_pgadb_relations(dbname)

    logger.info("Extract profiles from roary result matrix...")
    matrix_file = files.joinpath(output_dir, "roary", "gene_presence_absence.csv")
    profiles, total_isolates = extract_profiles(matrix_file, dbname)

    logger.info("Collecting allele profiles and making allele frequencies and reference sequence...")
    ffn_dir = files.joinpath(output_dir, "FFN")
    profile_file = files.joinpath(output_dir, "allele_profiles.tsv")
    profiles, freq = collect_allele_infos(profiles, ffn_dir)

    logger.info("Checking duplicating loci by self-blastp...")
    blastp_out_file, ref_length = reference_self_blastp(output_dir, freq)

    logger.info("Filter out high identity loci and drop loci which occurrence less than {}...".format(drop_by_occur))
    filtered_loci = filter_locus(blastp_out_file, ref_length, total_isolates, drop_by_occur)
    os.remove(blastp_out_file)

    logger.info("Updating and saving profiles...")
    freq = {l: freq[l] for l in filtered_loci}
    profiles = profiles[profiles.index.isin(filtered_loci)]
    profiles.to_csv(profile_file, sep="\t")

    logger.info("Saving allele sequences...")
    save_sequences(freq, dbname)

    logger.info("Making dynamic schemes...")
    make_schemes(freq, total_isolates)
    logger.info("Done!!")
    return dbname
