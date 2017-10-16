import json
import os
from concurrent.futures import ProcessPoolExecutor
import functional
import pandas as pd
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline

from src.models import logs
from src.utils import files, seq, cmds, operations
from src.utils.db import load_database_config, sql_query

BLAST_COLUMNS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                 "sstart", "send", "evalue", "bitscore"]


def identify_loci(args):
    filename, out_dir = args
    os.system(cmds.form_prodigal_cmd(filename, out_dir))
    return filename


def profile_by_query(filename, genome_id, selected_loci):
    # TODO: collect new alleles from here
    allele_ids = [operations.make_seqid(rec.seq) for rec in SeqIO.parse(filename, "fasta")]
    allele_ids = ",".join("'{}'".format(x) for x in allele_ids)
    locus_ids = ",".join("'{}'".format(x) for x in selected_loci)
    query = "select locus_id, allele_id from sequence where allele_id in ({}) and locus_id in ({});".format(allele_ids, locus_ids)
    profile = sql_query(query).drop_duplicates("allele_id")  # ensure allele_id is mapped only once
    profile = profile.drop_duplicates("locus_id").set_index("locus_id")  # ensure locus_id exists only once
    profile = profile.rename(columns={"allele_id": genome_id}).iloc[:, 0]
    return profile


def identify_and_profile(args):
    filename, query_dir, selected_loci = args
    in_file = os.path.join(query_dir, filename)
    out_dir = os.path.join(query_dir, "temp")
    genome_id = files.fasta_filename(filename)
    target_file = os.path.join(out_dir, genome_id + ".locus.fna")
    identify_loci(in_file, out_dir)
    profile = profile_by_query(target_file, genome_id)
    profile = profile[profile.index.isin(selected_loci)]
    return profile


def rename(query_dir, input_dir):
    namemap = {}
    for i, filename in enumerate(sorted(os.listdir(input_dir)), 1):
        file = SeqIO.parse(files.joinpath(input_dir, filename), "fasta")
        records = []
        for j, record in enumerate(file, 1):
            newid = "Genome_{i}::Contig_{j}".format(**locals())
            records.append(seq.new_record(newid, str(record.seq)))

        newname = "Genome_{i}.fa".format(**locals())
        SeqIO.write(records, files.joinpath(query_dir, newname), "fasta")
        namemap[files.replace_ext(newname)] = files.replace_ext(filename)
    return namemap


def profile_loci(refseq_fna, query_dir, output_dir, aligcov_cut, identity, threads, flat_file=True):
    if flat_file:
        refseqlen = (functional.seq(SeqIO.parse(refseq_fna, "fasta"))
                     .map(lambda rec: (rec.id, len(rec.seq)))
                     .to_dict())
    else:
        query = """SELECT locus_id, length(sequence.dna_seq) as ref_len
                   FROM sequence
                   INNER JOIN scheme
                   ON sequence.locus_id=scheme.locus_id
                   AND sequence.allele_id=scheme.ref_id;"""
        t = sql_query(query)
        refseqlen = {row["locus_id"]: row["ref_len"] for i, row in t.iterrows()}

    args = [(x, query_dir, refseq_fna, refseqlen, aligcov_cut, identity)
            for x in os.listdir(query_dir)]
    with ProcessPoolExecutor(threads) as executor:
        collect = {k: v for k, v in executor.map(extract_locus, args)}

    refseqs = list(refseqlen.keys())
    series = []
    for cid, loci in collect.items():
        xs = [s in loci for s in refseqs]
        ser = pd.Series(xs, name=cid, index=refseqs)
        series.append(ser)
    table = pd.concat(series, axis=1).sort_index(axis=0).sort_index(axis=1)
    table.to_csv(files.joinpath(output_dir, "locus_profiles.tsv"), sep="\t")


def extract_locus(args):
    filename, query_dir, refseq_fna, refseqlen, aligcov_cut, identity = args
    contig_file = files.joinpath(query_dir, filename)
    contig_id = filename.split(".")[0]

    db_dir = os.path.join(query_dir, contig_id)
    blastn_out_file = files.joinpath(query_dir, "{}.out".format(contig_id))

    compile_blastdb(contig_file, db_dir)
    query_db(refseq_fna, db_dir, blastn_out_file, BLAST_COLUMNS)
    matched_loci = identify_locus(blastn_out_file, refseqlen, aligcov_cut, identity, BLAST_COLUMNS)

    os.remove(db_dir + ".nhr")
    os.remove(db_dir + ".nin")
    os.remove(db_dir + ".nsq")
    os.remove(blastn_out_file)
    return contig_id, matched_loci


def compile_blastdb(input_file, output_file):
    cmd = "makeblastdb -in {} -dbtype nucl -out {}".format(input_file, output_file)
    os.system(cmd)


def query_db(query, db_dir, output_file, cols, threads=2):
    NcbiblastnCommandline(query=query, db=db_dir, out=output_file,
                          outfmt="'6 {}'".format(" ".join(cols)), num_threads=threads)()


def identify_locus(blast_out, seqlen, aligcov_cut, identity, cols):
    result = pd.read_csv(blast_out, sep="\t", header=None, names=cols)
    result["qlen"] = [seqlen[x] for x in result["qseqid"]]
    result["aligcov"] = (result["length"] - result["gapopen"]) / result["qlen"]
    result = result[(result["aligcov"] >= aligcov_cut) & (result["pident"] >= identity)]
    return set(result["qseqid"])


def exactly_match_in(records1, records2):
    for r1 in records1:
        r1_rev = r1.seq.reverse_complement()
        for r2 in records2:
            if r1.seq in r2.seq:
                return r1.id
            if r1_rev in r2.seq:
                return r1.id
    return None


def profile_alleles(query_dir, db_dir, output_dir, threads, occr_level, selector=None):
    locusfiles = files.joinpath(db_dir, "locusfiles")
    profiles = select_loci(db_dir, output_dir, occr_level, selector)

    collect = []
    with ProcessPoolExecutor(threads) as executor:
        for contig, profile in profiles.iteritems():
            contig_file = files.joinpath(query_dir, "{}.fa".format(contig))
            records = list(SeqIO.parse(contig_file, "fasta"))
            matched = profile[profile]

            args = [(locus, locusfiles, records) for locus in matched.index]
            series = pd.Series(name=contig)
            for x in executor.map(match_allele, args):
                if x:
                    locus, allele = x
                    series = series.set_value(locus, allele)
            collect.append(series)
    result = pd.concat(collect, axis=1)
    result.to_csv(files.joinpath(output_dir, "wgmlst.tsv"), sep="\t")


def select_loci(db_dir, output_dir, occr_level, selector):
    profile_file = files.joinpath(output_dir, "locus_profiles.tsv")
    profiles = pd.read_csv(profile_file, sep="\t", index_col=0)
    scheme = pd.read_csv(files.joinpath(db_dir, "scheme.tsv"), usecols=[0, 1], sep="\t")
    if not selector:
        selected_loci = scheme[scheme["occurence"] >= occr_level]["locus"]
    elif type(selector) == list:
        selected_loci = selector
    else:
        selected_loci = scheme["locus"]
    profiles = profiles[profiles.index.isin(selected_loci)]
    return profiles


def match_allele(args):
    locus, locusfiles, records = args
    alleles_file = files.joinpath(locusfiles, "{}.fa".format(locus))
    alleles = list(SeqIO.parse(alleles_file, "fasta"))
    matched_allele = exactly_match_in(alleles, records)
    if matched_allele:
        return locus, matched_allele
    return None


def profiling(output_dir, input_dir, db_dir, occr_level, threads,
              logger=None, aligcov_cut=0.5, identity=90, flat_file=True):
    load_database_config()
    if not logger:
        logger = logs.console_logger(__name__)

    logger.info("Renaming contigs...")
    query_dir = files.joinpath(output_dir, "query")
    files.create_if_not_exist(query_dir)
    namemap = rename(query_dir, input_dir)
    with open(files.joinpath(output_dir, "namemap.json"), "w") as f:
        f.write(json.dumps(namemap))

    if flat_file:
        logger.info("Profiling loci...")
        refseq_fna = files.joinpath(db_dir, "panRefSeq.fa")
        profile_loci(refseq_fna, query_dir, output_dir, aligcov_cut, identity, threads)

        logger.info("Allocating alleles...")
        profile_alleles(query_dir, db_dir, output_dir, threads, occr_level)
    else:
        logger.info("Identifying loci and allocating alleles...")

        # select loci by scheme
        query = "select locus_id from scheme where occurence>={};".format(occr_level)
        selected_loci = set(sql_query(query).iloc[:, 0])

        temp_dir = os.path.join(query_dir, "temp")
        files.create_if_not_exist(temp_dir)

        collect = []
        args = [(os.path.join(query_dir, filename), temp_dir) for filename in os.listdir(query_dir) if filename.endswith(".fa")]
        with ProcessPoolExecutor(threads) as executor:
            for filename in executor.map(identify_loci, args):
                genome_id = files.fasta_filename(filename)
                target_file = os.path.join(temp_dir, genome_id + ".locus.fna")
                profile = profile_by_query(target_file, genome_id, selected_loci)
                collect.append(profile)
        result = pd.concat(collect, axis=1)
        result.to_csv(files.joinpath(output_dir, "wgmlst.tsv"), sep="\t")
