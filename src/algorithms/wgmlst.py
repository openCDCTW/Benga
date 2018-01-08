import json
import os
import shutil
import functional
import pandas as pd
import functools
from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO
import subprocess
from src.models import logs
from src.utils import files, seq, cmds, operations
from src.utils.db import load_database_config, from_sql, append_to_sql

MLST = ["aroC_1", "aroC_2", "aroC_3", "dnaN", "hemD", "hisD", "purE", "sucA_1", "sucA_2", "thrA_2", "thrA_3"]
virulence_genes = ["lpfA", "lpfA_1", "lpfA_2", "lpfA_3", "lpfA_4", "lpfB", "lpfB_1", "lpfB_2", "lpfC", "lpfC_1",
                   "lpfD", "lpfD_1", "lpfD_2", "lpfE", # "lpfC''", "lpfC''_1", "lpfC''_2", "lpfC'_3",
                   "fimA_1", "fimA_2", "fimA_4", "fimA_5", "fimA_6", "fimC_1", "fimC_2", "fimC_3",
                   "fimD_1", "fimD_2", "fimD_3", "fimD_4", "fimD_5", "fim_2", "viaA_1", "viaA_2",
                   "fur_1", "fur_2", "rpoS", "rpoS_1", "rpoS_2", "spvB", "spvB_1", "spvB_2", "spvC"]


def identify_loci(args):
    filename, out_dir = args
    subprocess.run(cmds.form_prodigal_cmd(filename, out_dir), shell=True)
    genome_id = files.fasta_filename(filename)
    target_file = os.path.join(out_dir, genome_id + ".locus.fna")
    alleles = {operations.make_seqid(rec.seq): (rec.seq, rec.seq.translate(table=11))
               for rec in SeqIO.parse(target_file, "fasta")}
    return genome_id, alleles


def profile_by_query(alleles, genome_id, selected_loci, database):
    locus_ids = ",".join("'{}'".format(x) for x in selected_loci)
    allele_ids = ",".join("'{}'".format(x) for x in alleles.keys())
    query = "select allele_id, locus_id " \
            "from alleles " \
            "where allele_id in ({}) and locus_id in ({});".format(allele_ids, locus_ids)
    # ensure allele_id is mapped only once
    profile = from_sql(query, database=database).drop_duplicates("allele_id")
    # ensure locus_id exists only once
    profile = profile.drop_duplicates("locus_id").set_index("locus_id")
    profile = profile.rename(columns={"allele_id": genome_id}).iloc[:, 0]
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


def profile_loci(refseq_fna, query_dir, output_dir, aligcov_cut, identity, threads):
    refseqlen = (functional.seq(SeqIO.parse(refseq_fna, "fasta"))
                 .map(lambda rec: (rec.id, len(rec.seq)))
                 .to_dict())

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

    seq.compile_blastndb(contig_file, db_dir)
    seq.query_blastndb(refseq_fna, db_dir, blastn_out_file, seq.BLAST_COLUMNS)
    matched_loci = identify_locus(blastn_out_file, refseqlen, aligcov_cut, identity, seq.BLAST_COLUMNS)

    os.remove(db_dir + ".nhr")
    os.remove(db_dir + ".nin")
    os.remove(db_dir + ".nsq")
    os.remove(blastn_out_file)
    return contig_id, matched_loci


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
        selected_loci = scheme[scheme["occurrence"] >= occr_level]["locus"]
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


def make_ref_blastpdb(ref_db_file, database):
    query = "select loci.locus_id, alleles.peptide_seq " \
            "from loci inner join alleles " \
            "on loci.ref_allele = alleles.allele_id;"
    refs = from_sql(query, database=database)

    ref_recs = [seq.new_record(row["locus_id"], row["peptide_seq"], seqtype="protein") for _, row in refs.iterrows()]
    ref_fasta = ref_db_file + ".fasta"
    seq.save_records(ref_recs, ref_fasta)

    seq.compile_blastpdb(ref_fasta, ref_db_file)
    os.remove(ref_fasta)


def blast_for_new_alleles(candidates, alleles, ref_db, temp_dir, identity=95):
    '''blastp for locus with 95% identity and coverage 90%'''
    filename = "new_allele_candidates"
    candidate_file = os.path.join(temp_dir, filename + ".fasta")
    recs = [seq.new_record(cand, alleles[cand][1], seqtype="protein") for cand in candidates]
    seq.save_records(recs, candidate_file)

    blastp_out_file = files.joinpath(temp_dir, "{}.blastp.out".format(filename))
    seq.query_blastpdb(candidate_file, ref_db, blastp_out_file, seq.BLAST_COLUMNS, cov=90)

    blastp_out = pd.read_csv(blastp_out_file, sep="\t", header=None, names=seq.BLAST_COLUMNS)
    blastp_out = blastp_out[blastp_out["pident"] >= identity].drop_duplicates("qseqid")
    new_alleles = {row["qseqid"]: row["sseqid"] for _, row in blastp_out.iterrows()}
    return new_alleles


def update_database(new_allels, alleles):
    cols = ["allele_id", "locus_id", "dna_seq", "peptide_seq", "count"]
    collect = []
    for allele_id, locus_id in new_allels.items():
        dna = str(alleles[allele_id][0])
        peptide = str(alleles[allele_id][1])
        count = 1  # TODO: not real counting number
        collect.append((allele_id, locus_id, dna, peptide, count))
    result = pd.DataFrame(collect, columns=cols)
    append_to_sql("alleles", result)
    return result[["allele_id", "locus_id"]]


def add_new_alleles(id_allele_list, ref_db, temp_dir):
    all_alleles = functools.reduce(lambda x, y: {**x, **y[1]}, id_allele_list, {})
    existed_alleles = from_sql("select allele_id from alleles;")
    candidates = list(filter(lambda x: x not in existed_alleles["allele_id"], all_alleles.keys()))
    new_alleles = blast_for_new_alleles(candidates, all_alleles, ref_db, temp_dir)
    if new_alleles.keys():
        update_database(new_alleles, all_alleles)


def profiling(output_dir, input_dir, database, threads, occr_level=None, selected_loci=None, logger=None):
    load_database_config()
    if not logger:
        logger = logs.console_logger(__name__)

    logger.info("Renaming contigs...")
    query_dir = files.joinpath(output_dir, "query")
    files.create_if_not_exist(query_dir)
    namemap = rename(query_dir, input_dir)
    with open(files.joinpath(output_dir, "namemap.json"), "w") as f:
        f.write(json.dumps(namemap))

    logger.info("Identifying loci and allocating alleles...")
    if selected_loci:
        selected_loci = set(selected_loci)
    else:  # select loci by scheme
        query = "select locus_id from loci where occurrence>={};".format(occr_level)
        selected_loci = set(from_sql(query, database=database).iloc[:, 0])

    temp_dir = os.path.join(query_dir, "temp")
    files.create_if_not_exist(temp_dir)
    ref_db = os.path.join(temp_dir, "ref_blastpdb")
    make_ref_blastpdb(ref_db, database)

    args = [(os.path.join(query_dir, filename), temp_dir)
            for filename in os.listdir(query_dir) if filename.endswith(".fa")]
    with ProcessPoolExecutor(threads) as executor:
        id_allele_list = list(executor.map(identify_loci, args))
    add_new_alleles(id_allele_list, ref_db, temp_dir)

    collect = []
    for genome_id, alleles in id_allele_list:
        profile = profile_by_query(alleles, genome_id, selected_loci, database)
        collect.append(profile)
    result = pd.concat(collect, axis=1)
    result.to_csv(files.joinpath(output_dir, "wgmlst.tsv"), sep="\t")

    shutil.rmtree(query_dir)


def mlst_profiling(output_dir, input_dir, database, threads, logger=None):
    return profiling(output_dir, input_dir, database, threads, logger=logger, selected_loci=MLST)


def virulence_profiling(output_dir, input_dir, database, threads, logger=None):
    return profiling(output_dir, input_dir, database, threads, logger=logger, selected_loci=virulence_genes)
