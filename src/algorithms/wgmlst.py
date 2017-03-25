import json
import os
from collections import defaultdict
import functional
import pandas as pd
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from concurrent.futures import ProcessPoolExecutor

from src.utils import files, seq
from src.models import logs


BLAST_COLUMNS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                 "sstart", "send", "evalue", "bitscore"]


def rename(assemble_dir, query_dir):
    namemap = {}
    for i, filename in enumerate(sorted(os.listdir(query_dir)), 1):
        file = SeqIO.parse(files.joinpath(query_dir, filename), "fasta")
        records = []
        for j, record in enumerate(file, 1):
            newid = "Assembly_{i}::Contig_{j}".format(**locals())
            records.append(seq.new_record(newid, str(record.seq)))

        newname = "Assembly_{i}.fa".format(**locals())
        SeqIO.write(records, files.joinpath(assemble_dir, newname), "fasta")
        namemap[newname.split(".")[0]] = filename.split(".")[0]
    return namemap


def profile_loci(refseq_fna, assemble_dir, output_dir, refseqdb_dir, aligcov_cut, identity, threads):
    refseqlen = (functional.seq(SeqIO.parse(refseq_fna, "fasta"))
                 .map(lambda rec: (rec.id, len(rec.seq)))
                 .to_dict())

    # TODO: needs refactor
    args = [(assemble_dir, output_dir, refseqdb_dir, aligcov_cut, identity, refseqlen, x)
            for x in os.listdir(assemble_dir)]

    collect = {}
    compile_blastdb(refseq_fna, refseqdb_dir)
    with ProcessPoolExecutor(threads) as executor:
        for k, v in executor.map(extract_locus, args):
            collect[k] = v

    refseqs = list(refseqlen.keys())
    series = []
    for cid, loci in collect.items():
        xs = [s in loci for s in refseqs]
        ser = pd.Series(xs, name=cid, index=refseqs)
        series.append(ser)
    table = pd.concat(series, axis=1).sort_index(axis=0).sort_index(axis=1)
    table.to_csv(files.joinpath(output_dir, "locusAP.tsv"), sep="\t")


def extract_locus(args):
    assemble_dir, output_dir, refseqdb_dir, aligcov_cut, identity, refseqlen, filename = args
    contig_file = files.joinpath(assemble_dir, filename)
    contig_id = filename.split(".")[0]
    blastn_out_file = files.joinpath(output_dir, "{}.out".format(contig_id))
    query_db(contig_file, refseqdb_dir, blastn_out_file, BLAST_COLUMNS)
    matched_loci = identify_locus(blastn_out_file, refseqlen, aligcov_cut, identity, BLAST_COLUMNS)
    os.remove(blastn_out_file)
    return contig_id, matched_loci


def compile_blastdb(input_file, output_dir):
    cmd = "makeblastdb -in {} -dbtype nucl -out {}".format(input_file, output_dir)
    os.system(cmd)


def query_db(query, db_dir, output_file, cols, threads=2):
    NcbiblastnCommandline(query=query, db=db_dir, out=output_file,
                          outfmt="'6 {}'".format(" ".join(cols)), num_threads=threads)()


def identify_locus(blast_out, seqlen, aligcov_cut, identity, cols):
    result = pd.read_csv(blast_out, sep="\t", header=None, names=cols)
    result["slen"] = [seqlen[x] for x in result["sseqid"]]
    result["aligcov"] = (result["length"] - result["gapopen"]) / result["slen"]
    result = result[(result["aligcov"] >= aligcov_cut) & (result["pident"] >= identity)]
    return set(result["sseqid"])


def maxlen_locus(locus, pair):
    if all(map(lambda x: x == 0, pair)):
        return locus, 0
    else:
        # only access the last of max length allele number in values
        return locus, max(reversed(pair), key=lambda x: x[1])[0]


def exactly_match_in(records1, records2):
    for r1 in records1:
        for r2 in records2:
            if r1.seq in r2.seq:
                return r1.id
    return None


def profile_alleles(assemble_dir, db_dir, output_dir, threads, occr_level=90, selector=None):
    locusfiles = files.joinpath(db_dir, "locusfiles")
    profile_file = files.joinpath(output_dir, "locusAP.tsv")

    # select loci to profile depends on scheme
    scheme = pd.read_csv(files.joinpath(db_dir, "scheme.csv"), usecols=[0, 1])
    profiles = pd.read_csv(profile_file, sep="\t", index_col=0)
    if not selector:
        selected_loci = scheme[scheme["occurence"] >= occr_level]["locus"]
    elif type(selector) == list:
        selected_loci = selector
    else:
        selected_loci = scheme["locus"]
    profiles = profiles[profiles.index.isin(selected_loci)]

    # profiling
    collect = []
    with ProcessPoolExecutor(threads) as executor:
        for contig, profile in profiles.iteritems():
            contig_file = files.joinpath(assemble_dir, "{}.fa".format(contig))
            records = list(SeqIO.parse(contig_file, "fasta"))
            matched = profile[profile]

            args = [(locus, locusfiles, records) for locus in matched.index]
            series = pd.Series(name=contig)
            for x in executor.map(run, args):
                if x:
                    locus, allele = x
                    series = series.set_value(locus, allele)
            collect.append(series)
    result = pd.concat(collect, axis=1)
    result.to_csv(files.joinpath(output_dir, "wgmlst.tsv"), sep="\t", index=False)


def run(args):
    locus, locusfiles, records = args
    alleles_file = files.joinpath(locusfiles, "{}.fa".format(locus))
    alleles = list(SeqIO.parse(alleles_file, "fasta"))
    matched_allele = exactly_match_in(alleles, records)
    if matched_allele:
        return locus, matched_allele
    return None


def allocate_alleles(assemble_dir, db_dir, output_dir):
    profile_file = files.joinpath(output_dir, "locusAP.tsv")
    profiles = pd.read_csv(profile_file, sep="\t", index_col=0)

    for qfn, profile in profiles.iteritems():



    for filename in os.listdir(assemble_dir):
        qfn = filename.split(".")[0]
        os.mkdir(files.joinpath(output_dir, qfn))

        blastnfile = files.joinpath(assemble_dir, filename)
        compile_blastdb(blastnfile, files.joinpath(output_dir, "AssemblyDB_{}".format(qfn)))

        allprofile = defaultdict(list)
        newloc = []
        aldic = {}
        # TODO: parallelize
        # parse each contig profile
        for line in open(files.joinpath(output_dir, "locusAP.{}.list".format(qfn)), "r").read().splitlines():
            locus, hits = line.split("\t")
            if hits == "0":
                allprofile[locus] = [0]  # set locus hit to 0
                newloc.append(locus)  # potential new locus
            else:
                # parse refseq locus for dict from allele id to allele length
                qfile = files.joinpath(db_dir, "locusfiles", "{}.fa.new".format(locus))
                for record in SeqIO.parse(qfile, "fasta"):
                    alleleno = record.id.split("::")[1]
                    aldic[(locus, alleleno)] = len(record.seq)
                # append refseq locus for query
                os.system("cat {} >> {}/qfile.fa".format(qfile, output_dir))

        # query contig blastdb using collected refseq locus
        query_db(query=files.joinpath(output_dir, "qfile.fa"),
                 db_dir=files.joinpath(output_dir, "AssemblyDB_{}".format(qfn)),
                 output_file=files.joinpath(output_dir, "blast.{}.out".format(qfn)),
                 cols=BLAST_COLUMNS)
        os.remove(files.joinpath(output_dir, "qfile.fa"))

        # parse blast output
        blast_out = files.joinpath(output_dir, "blast.{}.out".format(qfn))
        blast = pd.read_csv(blast_out, sep="\t", header=None, names=BLAST_COLUMNS)
        blast["loc"] = blast["qseqid"].str.split("::").str[0]
        blast["loc"] = blast["loc"].str.split(".").str[0]
        blast["allno"] = blast["qseqid"].str.split("::").str[1]

        # set criteria
        blast = blast[(blast["pident"] == 100) & (blast["mismatch"] == 0) & (blast["gapopen"] == 0)]
        blast = blast[blast["length"] == aldic[(blast["loc"], blast["allno"])]]
        blast = blast[["loc", "allno", "length"]]
        for locus, df in blast.groupby("loc"):
            allprofile[locus].append(list(df.to_records(index=False)))

        # output profile
        output_file = files.joinpath(output_dir, "allele.{}.profile".format(qfn))
        (functional.seq(allprofile.items())
         .map(lambda locus, pair: maxlen_locus(locus, pair))
         .sorted(key=lambda x: x[0])
         .map(lambda x: x[0] + "\t" + str(x[1]))
         .to_file(output_file, delimiter="\n"))


def make_scheme_profile(db_dir, output_dir, occr_level=90):
    scheme = pd.read_csv(files.joinpath(db_dir, "scheme.csv"), usecols=[0, 1])
    selected_loci = scheme[scheme["occurence"] >= occr_level]["locus"]

    # select loci in profiles
    # output profiles

    for sch, loci in scheme.items():
        sink = files.joinpath(output_dir, "wgMLST_{}.tsv".format(sch))
        # table.loc[loci, :].to_csv(sink, sep="\t")


def profiling(output_dir, input_dir, db_dir, threads, logger=None, aligcov_cut=0.5, identity=90):
    if not logger:
        logger = logs.console_logger(__name__)

    logger.info("Renaming contigs...")
    assemble_dir = files.joinpath(output_dir, "query_assembly")
    files.create_if_not_exist(assemble_dir)
    namemap = rename(assemble_dir, input_dir)
    with open(files.joinpath(output_dir, "namemap.json"), "w") as f:
        f.write(json.dumps(namemap))

    logger.info("Profiling loci...")
    refseq_fna = files.joinpath(db_dir, "panRefSeq.fa")
    refseqdb_dir = files.joinpath(db_dir, "panRefSeq")
    profile_loci(refseq_fna, assemble_dir, output_dir, refseqdb_dir, aligcov_cut, identity, threads)

    logger.info("Allocating alleles...")
    profile_alleles(assemble_dir, db_dir, output_dir, threads)

    logger.info("Collecting profiles by scheme...")
    # make_scheme_profile(assemble_dir, db_dir, output_dir)

    logger.info("Output...")
    # wgmlst_dir = files.joinpath(scheme_dir, "wgMLST")
    # files.clear_folder(wgmlst_dir)
    # os.system("mv {}/wgMLST_* {}".format(output_dir, wgmlst_dir))
    # os.system("rm -rf Assembly*")

