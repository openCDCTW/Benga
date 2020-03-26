import functools
import os
import shutil
import uuid
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
import sqlalchemy as sa
from Bio import SeqIO

from src.algorithms.bionumerics import to_bionumerics_format
from src.utils import files, cmds, operations, logs, seq
from src.utils import db
from src.utils.alleles import filter_duplicates

CURRENT_SCRIPT_DIRECTORY = os.path.dirname(os.path.realpath(__file__))
MODELS_DIRECTORY = os.path.abspath(os.path.join(CURRENT_SCRIPT_DIRECTORY, "..", "..", 'models'))


def identify_alleles(genome_file, outdir, model):
    """Identify CDS from the outcome of prodigal."""
    genome_id = files.get_fileroot(genome_file)
    training_file = os.path.join(MODELS_DIRECTORY, model + '.trn')
    outfile = os.path.join(outdir, genome_id + '.locus.fna')
    prodigal_cmd = cmds.form_prodigal_cmd(genome_file, outfile, training_file)
    cmds.execute_cmd(prodigal_cmd)
    alleles = OrderedDict()
    for rec in SeqIO.parse(outfile, "fasta"):
        allele_id = operations.make_seqid(rec.seq)
        content = (rec.seq, rec.seq.translate(table=11))
        alleles[allele_id] = content
    return genome_id, alleles


def profile_by_query(alleles, genome_id, selected_loci, database):
    """Profiling a genome by query database with allele id.
    Ensure an allele is mapped to a locus.
    """
    query = sa.select([db.PAIRS]).where(sa.and_(
        db.PAIRS.c['allele_id'].in_(alleles.keys()),
        db.PAIRS.c['locus_id'].in_(selected_loci)
    ))
    # ensure allele_id is mapped only once
    profile = db.from_sql(query, database=database).drop_duplicates("allele_id")
    # ensure locus_id exists only once
    profile = profile.sort_values("allele_id", kind='mergesort').drop_duplicates("locus_id").set_index("locus_id")
    profile = profile.rename(columns={"allele_id": genome_id}).iloc[:, 0]
    return profile


def make_ref_blastpdb(ref_db_file, database):
    query = sa.select([db.LOCI.c["locus_id"], db.ALLELES.c["peptide_seq"]]).select_from(
        db.LOCI.join(db.ALLELES, db.LOCI.c["ref_allele"] == db.ALLELES.c["allele_id"])
    )
    refs = db.from_sql(query, database=database)

    ref_recs = [seq.new_record(row["locus_id"], row["peptide_seq"], seqtype="protein") for _, row in refs.iterrows()]
    ref_fasta = ref_db_file + ".fasta"
    seq.save_records(ref_recs, ref_fasta)

    seq.compile_blastpdb(ref_fasta, ref_db_file)
    os.remove(ref_fasta)


def blast_for_new_alleles(candidates, alleles, ref_db, temp_dir, pid, threads):
    """Blast unmapped alleles with existing locus and identify potentially new alleles."""
    filename = "new_allele_candidates_" + pid
    candidate_file = os.path.join(temp_dir, filename + ".fasta")
    recs = [seq.new_record(cand, alleles[cand][1], seqtype="protein") for cand in candidates]
    seq.save_records(recs, candidate_file)

    blastp_out_file = os.path.join(temp_dir, "{}.blastp.out".format(filename))
    seq.query_blastpdb(candidate_file, ref_db, blastp_out_file, seq.BLAST_COLUMNS, threads)

    blastp_out = filter_duplicates(blastp_out_file, identity=95)
    blastp_out = blastp_out.drop_duplicates("qseqid")
    new_allele_pairs = [(row["qseqid"], row["sseqid"]) for _, row in blastp_out.iterrows()]
    return new_allele_pairs


def update_database(new_allele_pairs, alleles):
    """Update new alleles to database."""
    collect = []
    for allele_id, locus_id in new_allele_pairs:
        dna = str(alleles[allele_id][0])
        peptide = str(alleles[allele_id][1])
        count = 0
        collect.append((allele_id, dna, peptide, count))
    collect = pd.DataFrame(collect, columns=["allele_id", "dna_seq", "peptide_seq", "count"]).drop_duplicates()
    db.table_to_sql("alleles", collect)
    pairs = pd.DataFrame(new_allele_pairs, columns=["allele_id", "locus_id"]).drop_duplicates()
    db.table_to_sql("pairs", pairs)
    return pairs


def add_new_alleles(id_allele_list, ref_db, temp_dir, pid, threads):
    """Identify new alleles and add them to database."""
    all_alleles = functools.reduce(lambda x, y: {**x, **y[1]}, id_allele_list, {})
    query = sa.select([db.ALLELES.c["allele_id"]]).where(db.ALLELES.c["allele_id"].in_(all_alleles.keys()))
    existed_alleles = db.from_sql(query)["allele_id"].to_list()
    candidates = list(filter(lambda x: x not in existed_alleles, all_alleles.keys()))
    new_allele_pairs = blast_for_new_alleles(candidates, all_alleles, ref_db, temp_dir, pid, threads)
    if new_allele_pairs:
        update_database(new_allele_pairs, all_alleles)


def profiling(output_dir, input_dir, database, threads, occr_level=None, selected_loci=None,
              profile_file="profile", enable_adding_new_alleles=True, generate_profiles=True,
              generate_bn=True, logger=None, debug=False):
    pid = uuid.uuid4().hex[0:8]
    if not logger:
        lf = logs.LoggerFactory()
        lf.addConsoleHandler()
        lf.addFileHandler(os.path.join(output_dir, "profiling_" + pid + ".log"))
        logger = lf.create()
    db.load_database_config(logger=logger)

    logger.info("Formating contigs...")
    query_dir = os.path.join(output_dir, "query_{}".format(pid))
    os.makedirs(query_dir, exist_ok=True)
    contighandler = files.ContigHandler()
    contighandler.format(input_dir, query_dir)

    logger.info("Selecting loci by specified scheme {}%...".format(occr_level))
    if selected_loci:
        selected_loci = set(selected_loci)
    else:  # select loci by scheme
        query = "select locus_id from loci where occurrence>={};".format(occr_level)
        selected_loci = set(db.from_sql(query, database=database).iloc[:, 0])

    logger.info("Making reference blastdb for blastp...")
    temp_dir = os.path.join(query_dir, "temp_{}".format(pid))
    os.makedirs(temp_dir, exist_ok=True)
    ref_db = os.path.join(temp_dir, "ref_blastpdb_{}".format(pid))
    make_ref_blastpdb(ref_db, database)

    logger.info("Identifying loci and allocating alleles...")
    to_do = []
    with ProcessPoolExecutor(threads) as executor:
        for filename in os.listdir(query_dir):
            if contighandler.is_fasta(filename):
                future = executor.submit(identify_alleles, os.path.join(query_dir, filename), temp_dir, database)
                to_do.append(future)
    done_iter = as_completed(to_do)
    id_allele_list = [done.result() for done in done_iter]

    if enable_adding_new_alleles:
        logger.info("Adding new alleles to database...")
        add_new_alleles(id_allele_list, ref_db, temp_dir, pid, threads)

    logger.info("Collecting allele profiles of each genomes...")
    if generate_profiles:
        to_do = []
        with ProcessPoolExecutor(threads) as executor:
            for genome_id, alleles in id_allele_list:
                future = executor.submit(profile_by_query, alleles, genome_id, selected_loci, database)
                to_do.append(future)
        done_iter = as_completed(to_do)
        collect = [done.result() for done in done_iter]
        result = pd.concat(collect, axis=1, sort=False)
        result.to_csv(os.path.join(output_dir, profile_file + ".tsv"), sep="\t")
        if generate_bn:
            bio = to_bionumerics_format(result)
            bio.to_csv(os.path.join(output_dir, "bionumerics_{}.csv".format(pid)), index=False)
    else:
        logger.info("Not going to output profiles.")

    if not debug and os.path.exists(query_dir):
        shutil.rmtree(query_dir)
    logger.info("Done!")
