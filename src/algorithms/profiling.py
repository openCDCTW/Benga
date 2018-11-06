import functools
import os
import re
import shutil
import subprocess
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
from Bio import SeqIO

from src.algorithms.bionumerics import to_bionumerics_format
from src.utils import files, cmds, operations, logs, seq
from src.utils.db import load_database_config, from_sql, table_to_sql, to_sql
from src.utils.alleles import filter_duplicates

MLST = ["aroC_1", "aroC_2", "aroC_3", "dnaN", "hemD", "hisD", "purE", "sucA_1", "sucA_2", "thrA_2", "thrA_3"]
virulence_genes = ["lpfA", "lpfA_1", "lpfA_2", "lpfA_3", "lpfA_4", "lpfB", "lpfB_1", "lpfB_2", "lpfC", "lpfC_1",
                   "lpfD", "lpfD_1", "lpfD_2", "lpfE", # "lpfC''", "lpfC''_1", "lpfC''_2", "lpfC'_3",
                   "fimA_1", "fimA_2", "fimA_4", "fimA_5", "fimA_6", "fimC_1", "fimC_2", "fimC_3",
                   "fimD_1", "fimD_2", "fimD_3", "fimD_4", "fimD_5", "fim_2", "viaA_1", "viaA_2",
                   "fur_1", "fur_2", "rpoS", "rpoS_1", "rpoS_2", "spvB", "spvB_1", "spvB_2", "spvC"]

# class Profiler(metaclass=abc.ABCMeta):
#     def __init__(self):
#         raise NotImplementedError()
#
# class BatchProfiler(Profiler):
#     def __init__(self):
#         self.__profile = None
#
#     def run(self, ):
#
#
# class BatchSelector(metaclass=abc.ABCMeta):
#     def __init__(self):
#         raise NotImplementedError()
#
# class BatchSchemeSelector(BatchSelector):
#     def __init__(self):
#         pass
#
#     def
#
# class BatchGeneSelector(BatchSelector):


def identify_alleles(args):
    filename, out_dir, model = args
    subprocess.run(cmds.form_prodigal_cmd(filename, out_dir, model), shell=True)
    genome_id = files.fasta_filename(filename)
    target_file = os.path.join(out_dir, genome_id + ".locus.fna")
    alleles = {operations.make_seqid(rec.seq): (rec.seq, rec.seq.translate(table=11))
               for rec in SeqIO.parse(target_file, "fasta")}
    return genome_id, alleles


def update_allele_counts(counter, database):
    table_to_sql("batch_add_counts", counter, database=database, append=False)
    query = "update alleles " \
            "set count = alleles.count + ba.count " \
            "from batch_add_counts as ba " \
            "where alleles.allele_id=ba.allele_id;"
    to_sql(query, database=database)
    to_sql("drop table batch_add_counts;", database=database)


def profile_by_query(alleles, genome_id, selected_loci, database):
    locus_ids = ",".join("'{}'".format(x) for x in selected_loci)
    allele_ids = ",".join("'{}'".format(x) for x in alleles.keys())
    query = "select allele_id, locus_id " \
            "from pairs " \
            "where allele_id in ({}) and locus_id in ({});".format(allele_ids, locus_ids)
    # ensure allele_id is mapped only once
    profile = from_sql(query, database=database).drop_duplicates("allele_id")
    # ensure locus_id exists only once
    profile = profile.drop_duplicates("locus_id").set_index("locus_id")
    profile = profile.rename(columns={"allele_id": genome_id}).iloc[:, 0]
    return profile


def generate_allele_len(recs):
    return {rec.id: len(rec.seq) for rec in recs}


def make_ref_blastpdb(ref_db_file, database):
    query = "select loci.locus_id, alleles.peptide_seq " \
            "from loci inner join alleles " \
            "on loci.ref_allele = alleles.allele_id;"
    refs = from_sql(query, database=database)

    ref_recs = [seq.new_record(row["locus_id"], row["peptide_seq"], seqtype="protein") for _, row in refs.iterrows()]
    ref_fasta = ref_db_file + ".fasta"
    seq.save_records(ref_recs, ref_fasta)
    ref_len = generate_allele_len(ref_recs)

    seq.compile_blastpdb(ref_fasta, ref_db_file)
    os.remove(ref_fasta)
    return ref_len


def blast_for_new_alleles(candidates, alleles, ref_db, temp_dir, ref_len):
    filename = "new_allele_candidates"
    candidate_file = os.path.join(temp_dir, filename + ".fasta")
    recs = [seq.new_record(cand, alleles[cand][1], seqtype="protein") for cand in candidates]
    seq.save_records(recs, candidate_file)
    allele_len = generate_allele_len(recs)

    blastp_out_file = files.joinpath(temp_dir, "{}.blastp.out".format(filename))
    seq.query_blastpdb(candidate_file, ref_db, blastp_out_file, seq.BLAST_COLUMNS)

    blastp_out = filter_duplicates(blastp_out_file, allele_len, ref_len, identity=95)
    blastp_out = blastp_out.drop_duplicates("qseqid")
    new_allele_pairs = [(row["qseqid"], row["sseqid"]) for _, row in blastp_out.iterrows()]
    return new_allele_pairs


def update_database(new_allele_pairs, alleles):
    collect = []
    for allele_id, locus_id in new_allele_pairs:
        dna = str(alleles[allele_id][0])
        peptide = str(alleles[allele_id][1])
        count = 0
        collect.append((allele_id, dna, peptide, count))
    collect = pd.DataFrame(collect, columns=["allele_id", "dna_seq", "peptide_seq", "count"]).drop_duplicates()
    table_to_sql("alleles", collect)
    pairs = pd.DataFrame(new_allele_pairs, columns=["allele_id", "locus_id"]).drop_duplicates()
    table_to_sql("pairs", pairs)
    return pairs


def add_new_alleles(id_allele_list, ref_db, temp_dir, ref_len):
    all_alleles = functools.reduce(lambda x, y: {**x, **y[1]}, id_allele_list, {})
    existed_alleles = from_sql("select allele_id from alleles;")["allele_id"].tolist()
    candidates = list(filter(lambda x: x not in existed_alleles, all_alleles.keys()))
    new_allele_pairs = blast_for_new_alleles(candidates, all_alleles, ref_db, temp_dir, ref_len)
    if new_allele_pairs:
        update_database(new_allele_pairs, all_alleles)


def profiling(output_dir, input_dir, database, threads, occr_level=None, selected_loci=None,
              enable_adding_new_alleles=True, generate_profiles=True, logger=None, debug=False):
    if not logger:
        lf = logs.LoggerFactory()
        lf.addConsoleHandler()
        lf.addFileHandler(files.joinpath(output_dir, "profiling.log"))
        logger = lf.create()
    load_database_config(logger=logger)

    logger.info("Formating contigs...")
    query_dir = files.joinpath(output_dir, "query")
    files.create_if_not_exist(query_dir)
    contighandler = files.ContigHandler()
    contighandler.new_format(input_dir, query_dir, replace_ext=True)
    namemap = contighandler.namemap

    model = re.search('[a-zA-Z]+\w[a-zA-Z]+', database).group(0)
    logger.info("Used model: {}".format(model))

    logger.info("Selecting loci by specified scheme {}%...".format(occr_level))
    if selected_loci:
        selected_loci = set(selected_loci)
    else:  # select loci by scheme
        query = "select locus_id from loci where occurrence>={};".format(occr_level)
        selected_loci = set(from_sql(query, database=database).iloc[:, 0])

    logger.info("Making reference blastdb for blastp...")
    temp_dir = os.path.join(query_dir, "temp")
    files.create_if_not_exist(temp_dir)
    ref_db = os.path.join(temp_dir, "ref_blastpdb")
    ref_len = make_ref_blastpdb(ref_db, database)

    logger.info("Identifying loci and allocating alleles...")
    args = [(os.path.join(query_dir, filename), temp_dir, model)
            for filename in os.listdir(query_dir) if filename.endswith(".fa")]
    with ThreadPoolExecutor(threads) as executor:
        id_allele_list = list(executor.map(identify_alleles, args))

    if enable_adding_new_alleles:
        logger.info("Adding new alleles to database...")
        add_new_alleles(id_allele_list, ref_db, temp_dir, ref_len)

    logger.info("Collecting allele profiles of each genomes...")
    allele_counts = Counter()
    if generate_profiles:
        collect = []
        for genome_id, alleles in id_allele_list:
            profile = profile_by_query(alleles, genome_id, selected_loci, database)
            collect.append(profile)
            allele_counts.update(alleles.keys())
        result = pd.concat(collect, axis=1)
        result.columns = list(map(lambda x: namemap[x], result.columns))
        result.to_csv(files.joinpath(output_dir, "profile.tsv"), sep="\t")
        bio = to_bionumerics_format(result)
        bio.to_csv(os.path.join(output_dir, 'bionumerics.csv'), index=False)
    else:
        logger.info("Not going to output profiles.")
        for genome_id, alleles in id_allele_list:
            allele_counts.update(alleles.keys())

    allele_counts = pd.DataFrame(allele_counts, index=[0]).T\
        .reset_index().rename(columns={"index": "allele_id", 0: "count"})
    update_allele_counts(allele_counts, database)
    if not debug:
        shutil.rmtree(query_dir)
    logger.info("Done!")
