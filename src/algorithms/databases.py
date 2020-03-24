import os
import subprocess
from functools import reduce
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor

import pandas as pd
from Bio import SeqIO

from .profiling import profiling
from src.utils import seq, files, cmds, operations, db, logs
from src.utils.alleles import filter_duplicates


def move_file(annotate_dir, dest_dir, ext):
    for folder in os.listdir(annotate_dir):
        path = os.path.join(annotate_dir, folder)
        for file in os.listdir(path):
            current_file = os.path.join(path, file)
            if file.endswith(ext):
                os.rename(current_file, os.path.join(dest_dir, file))


def filter_tRNA(matrix):
    """Filter out the locus of which description contains 'tRNA-***'."""
    return matrix[~matrix['description'].str.contains(r"tRNA-\w+\(\w{3}\)", regex=True)]


def filter_rRNA(matrix):
    """Filter out the locus of which description contains 'ribosomal RNA'."""
    return matrix[~matrix['description'].str.contains("ribosomal RNA&subunit")]


def extract_profiles(roary_matrix_file, dbname, metadata_cols=13):
    """Extract profiles from Roary outcome and identify the column names."""
    matrix = pd.read_csv(roary_matrix_file, low_memory=False)
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
    db.table_to_sql("locus_meta", meta, dbname)


def collect_allele_info(profiles, ffn_dir):
    """Collect allele information from ffn files.
    Collect all found loci, including paralogs, and generate allele id from allele sequence.
    Allele frequency and locus frequency are calculated.
    """
    freq = defaultdict(Counter)
    clusters = {}
    for rows in profiles.iterrows():
        clusters[rows[0]] = reduce(lambda x, y: x + y, [i.split('\t') for i in rows[1].dropna()], [])

    seqs = {}
    for ffn_file in os.listdir(ffn_dir):
        for record in SeqIO.parse(os.path.join(ffn_dir, ffn_file), "fasta"):
            seqs[record.id] = record.seq

    for locus, alleles in clusters.items():
        freq[locus].update(list(map(lambda x: seqs[x], alleles)))
    return freq


def reference_self_blastp(output_dir, freq, threads):
    """Blast reference alleles of locus with themselves."""
    ref_recs = [seq.new_record(locus, counter.most_common(1)[0][0].translate(table=11)) for locus, counter in freq.items()]
    ref_faa = os.path.join(output_dir, "ref_seq.faa")
    seq.save_records(ref_recs, ref_faa)

    ref_db = os.path.join(output_dir, "ref_db")
    seq.compile_blastpdb(ref_faa, ref_db)

    blastp_out_file = os.path.join(output_dir, "ref_db.blastp.out")
    seq.query_blastpdb(ref_faa, ref_db, blastp_out_file, seq.BLAST_COLUMNS, threads)
    return blastp_out_file


def identify_pairs(df):
    """Remove duplicated locus pairs in qseqid and sseqid."""
    sseqids = df["sseqid"].tolist()
    pairs = []
    existed = set()
    for _, row in df.iterrows():
        if row["qseqid"] in sseqids and not row["qseqid"] in existed:
            existed.update(row["qseqid"])
            existed.update(row["sseqid"])
            pairs.append((row["qseqid"], row["sseqid"]))
    return pairs


def select_drop_loci(df):
    """Drop locus id started with 'group_' and end with '_\d'."""
    drop1 = df[df["locus_id"].str.contains("group_")]["locus_id"]
    drop2 = df[df["locus_id"].str.contains(r"_\d$")]["locus_id"]
    return set(drop1) | set(drop2)


def collect_high_occurrence_loci(pairs, total_isolates, drop_by_occur):
    """Collect locus. If duplicated, keep the high occurrence one."""
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


def filter_locus(blastp_out_file, total_isolates, drop_by_occur):
    blastp_out = filter_duplicates(blastp_out_file, identity=95)
    pairs = identify_pairs(blastp_out)
    filtered_loci = collect_high_occurrence_loci(pairs, total_isolates, drop_by_occur)
    return filtered_loci


def to_allele_table(data, dbname):
    """Save allele information to allele table in database."""
    df = pd.DataFrame(data, columns=["allele_id", "dna_seq", "peptide_seq", "count"])
    df = df.groupby("allele_id").agg({"dna_seq": "first", "peptide_seq": "first", "count": "sum"})
    df.reset_index(inplace=True)
    db.table_to_sql("alleles", df, dbname)


def to_pair_table(data, dbname):
    """Save allele-locus pair information to pair table in database."""
    df = pd.DataFrame(data, columns=["allele_id", "locus_id"])
    db.table_to_sql("pairs", df, dbname)


def save_sequences(freq, refseqs, dbname):
    """Collect sequence information and save into database."""
    alleles = []
    pairs = []
    for locus, counter in freq.items():
        allele = refseqs[locus]
        count = 0
        dna_seq = str(allele)
        pept_seq = str(allele.translate(table=11))
        allele_id = operations.make_seqid(dna_seq)
        alleles.append((allele_id, dna_seq, pept_seq, count))
        pairs.append((allele_id, locus))
    to_allele_table(alleles, dbname)
    to_pair_table(pairs, dbname)


def make_schemes(refseqs, total_isolates):
    """Decide scheme and reference allele of a locus."""
    schemes = db.from_sql("select locus_id, num_isolates from locus_meta;")
    schemes = schemes[schemes["locus_id"].isin(refseqs.keys())]
    schemes["occurrence"] = list(map(lambda x: round(x/total_isolates * 100, 2), schemes["num_isolates"]))
    schemes["ref_allele"] = list(map(lambda x: refseqs[x], schemes["locus_id"]))
    schemes = schemes[["locus_id", "occurrence", "ref_allele"]]
    db.table_to_sql("loci", schemes)


def update_schemes(directory, database, threads):
    indir, outdir = os.path.join(directory, 'Genomes'), os.path.join(directory, 'Profile')
    os.makedirs(outdir, exist_ok=True)
    profiling(outdir, indir, database, threads, 0)
    profile = pd.read_csv(os.path.join(outdir, 'profile.tsv'), sep='\t', index_col=0, low_memory=False)
    isolates = profile.shape[1]
    new_schemes = {rows[0]: rows[1].notna().sum() for rows in profile.iterrows()}

    locus_meta = db.from_sql("select * from locus_meta;", database=database)
    locus_meta = locus_meta[locus_meta['locus_id'].isin(list(new_schemes))]
    locus_meta['num_isolates'] = locus_meta['locus_id'].map(new_schemes)

    loci = db.from_sql("select * from loci;", database=database)
    loci["occurrence"] = (loci['locus_id'].map(new_schemes)/isolates*100).round(2)
    loci = loci[loci["occurrence"].notna()]

    alleles = db.from_sql("select * from alleles;", database=database)
    alleles = alleles[alleles['allele_id'].isin(loci['ref_allele'])]

    pairs = db.from_sql("select * from pairs;", database=database)
    pairs = pairs[pairs['allele_id'].isin(loci['ref_allele'])]
    db.dropdb(database)
    db.createdb(database)
    db.create_pgadb_relations(database)
    db.table_to_sql("locus_meta", locus_meta, database)
    db.table_to_sql("alleles", alleles, database)
    db.table_to_sql("pairs", pairs, database)
    db.table_to_sql("loci", loci, database)


def annotate_configs(input_dir, output_dir, logger=None, threads=8, training_file=None):
    if not logger:
        lf = logs.LoggerFactory()
        lf.addConsoleHandler()
        lf.addFileHandler(os.path.join(output_dir, "annotation.log"))
        logger = lf.create()

    logger.info("Formating contigs...")
    genome_dir = os.path.join(output_dir, "Genomes")
    os.makedirs(genome_dir, exist_ok=True)
    contighandler = files.ContigHandler()
    contighandler.format(input_dir, genome_dir)

    logger.info("Annotating...")
    annotate_dir = os.path.join(output_dir, "Annotated")
    os.makedirs(annotate_dir, exist_ok=True)

    prokka_cmd = []
    for filename in os.listdir(genome_dir):
        genome_file = os.path.join(genome_dir, filename)
        outdir = os.path.join(annotate_dir, files.get_fileroot(filename))
        prokka_cmd.append(cmds.form_prokka_cmd(genome_file=genome_file, outdir=outdir, training_file=training_file))
    with ProcessPoolExecutor(int(threads / 2)) as executor:
        executor.map(cmds.execute_cmd, prokka_cmd)

    logger.info("Moving protein CDS (.ffn) files...")
    ffn_dir = os.path.join(output_dir, "FFN")
    os.makedirs(ffn_dir, exist_ok=True)
    move_file(annotate_dir, ffn_dir, ".ffn")

    logger.info("Moving annotation (.gff) files...")
    gff_dir = os.path.join(output_dir, "GFF")
    os.makedirs(gff_dir, exist_ok=True)
    move_file(annotate_dir, gff_dir, ".gff")


def make_database(output_dir, drop_by_occur, logger=None, threads=2):
    if not logger:
        lf = logs.LoggerFactory()
        lf.addConsoleHandler()
        lf.addFileHandler(os.path.join(output_dir, "make_database.log"))
        logger = lf.create()
    db.load_database_config(logger=logger)

    logger.info("Calculating the pan genome...")
    min_identity = 95
    c = cmds.form_roary_cmd(os.path.join(output_dir, "GFF"), output_dir, min_identity, threads)
    logger.info("Run roary with following command: " + c)
    subprocess.run(c, shell=True)

    logger.info("Creating database")
    dbname = os.path.basename(output_dir[:-1] if output_dir.endswith("/") else output_dir)
    db.createdb(dbname)
    db.create_pgadb_relations(dbname)

    logger.info("Extract profiles from roary result matrix...")
    matrix_file = os.path.join(output_dir, "roary", "gene_presence_absence.csv")
    profiles, total_isolates = extract_profiles(matrix_file, dbname)

    logger.info("Collecting allele profiles and making allele frequencies and reference sequence...")
    ffn_dir = os.path.join(output_dir, "FFN")
    freq = collect_allele_info(profiles, ffn_dir)

    logger.info("Checking duplicated loci by self-blastp...")
    blastp_out_file = reference_self_blastp(output_dir, freq, threads)

    logger.info("Filter out high identity loci and drop loci which occurrence less than {}...".format(drop_by_occur))
    filtered_loci = filter_locus(blastp_out_file, total_isolates, drop_by_occur)
    os.remove(blastp_out_file)

    logger.info("Updating and saving profiles...")
    freq = {l: freq[l] for l in filtered_loci}

    logger.info("Saving allele sequences...")
    refseqs = {locus: counter.most_common(1)[0][0] for locus, counter in freq.items()}
    save_sequences(freq, refseqs, dbname)

    logger.info("Making dynamic schemes...")
    refseqs = dict(map(lambda x: (x[0], operations.make_seqid(x[1])), refseqs.items()))
    make_schemes(refseqs, total_isolates)

    logger.info("Update schemes...")
    update_schemes(output_dir, dbname, threads)
    logger.info("Done!!")
    return dbname
