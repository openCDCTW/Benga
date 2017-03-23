import json
import os
from collections import defaultdict
import functional
import pandas as pd
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline

from src.utils import files, seq
from src.models import logs


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


def profile_loci(refseq_fna, assemble_dir, output_dir, refseqdb_dir, id_locus, blast_cols):
    refseqlen = (functional.seq(SeqIO.parse(refseq_fna, "fasta"))
                 .map(lambda rec: (rec.id, len(rec.seq)))
                 .to_dict())

    collect = {}
    compile_blastdb(refseq_fna, refseqdb_dir)
    for filename in os.listdir(assemble_dir):
        contig_file = files.joinpath(assemble_dir, filename)
        contig_id = filename.split(".")[0]
        blastn_out = files.joinpath(output_dir, contig_id + ".out")

        query_db(contig_file, refseqdb_dir, blastn_out, blast_cols)
        collect[contig_id] = id_locus(blastn_out, refseqlen)

        os.remove(blastn_out)

    refseqs = list(refseqlen.keys())
    series = []
    for cid, loci in collect.items():
        xs = [s in loci for s in refseqs]
        ser = pd.Series(xs, name=cid, index=refseqs)
        series.append(ser)
    table = pd.concat(series, axis=1).sort_index(axis=0)
    table.to_csv(files.joinpath(output_dir, "locusAP.tsv"), sep="\t")


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


def allocate_alleles(assemble_dir, db_dir, output_dir, bn, blast_cols):
    for filename in os.listdir(assemble_dir):
        qfn = filename.split(".")[0]
        os.mkdir(files.joinpath(output_dir, qfn))

        blastnfile = files.joinpath(assemble_dir, filename)
        compile_blastdb(blastnfile, files.joinpath(output_dir, "AssemblyDB_" + qfn))

        allprofile = defaultdict(list)
        newloc = []
        aldic = {}
        # TODO: parallelize
        for line in open(files.joinpath(output_dir, "locusAP." + qfn + ".list"), "r").read().splitlines():
            locus, hits = line.split("\t")
            if hits == "0":
                allprofile[locus] = [0]
                newloc.append(locus)
            else:
                qfile = files.joinpath(db_dir, "locusfiles", locus + ".fa.new")
                for record in SeqIO.parse(qfile, "fasta"):
                    alleleno = record.id.split("::")[1]
                    aldic[(locus, alleleno)] = len(record.seq)
                os.system("cat {qfile} >> {output_dir}/qfile.fa".format(**locals()))

        os.system(bn(files.joinpath(output_dir, "qfile.fa"), files.joinpath(output_dir, "AssemblyDB_" + qfn),
                     files.joinpath(output_dir, "blast." + qfn + ".out")))
        os.remove(files.joinpath(output_dir, "qfile.fa"))
        blast = pd.read_csv(files.joinpath(output_dir, "blast." + qfn + ".out"), sep="\t", header=None, names=blast_cols)
        blast["loc"] = blast["qseqid"].str.split("::").str[0]
        blast["loc"] = blast["loc"].str.split(".").str[0]
        blast["allno"] = blast["qseqid"].str.split("::").str[1]

        for _, row in blast[(blast["pident"] == 100) & (blast["mismatch"] == 0) & (blast["gapopen"] == 0)].iterrows():
            if row["length"] == aldic[(row["loc"], row["allno"])]:
                allprofile[row["loc"]].append((row["allno"], row["length"]))

        # Find new alleles
        new_alleles = []
        records = []
        assdic = {record.id: str(record.seq) for record in SeqIO.parse(blastnfile, "fasta")}
        for _, row in blast[~blast["loc"].isin(allprofile.keys())].iterrows():
            allprofile[row["loc"]] = [0]
            new_alleles.append(row["loc"])
            records.append(seq.new_record(row["loc"], assdic[row["sseqid"]]))

        SeqIO.write(files.drop_duplicate(records, lambda x: x.id), files.joinpath(output_dir, qfn, "tot.1.new"), "fasta")
        with open(files.joinpath(output_dir, qfn, "tot.new"), "w") as file:
            file.write("\n".join(files.drop_duplicate(new_alleles)))

        candidate_loci = [maxlen_locus(locus, pair) for locus, pair in allprofile.items()]
        (functional.seq(sorted(candidate_loci, key=lambda x: x[0]))
         .map(lambda x: x[0] + "\t" + str(x[1]))
         .to_file(files.joinpath(output_dir, "allele." + qfn + ".profile"), delimiter="\n"))


def make_scheme_profile(assemble_dir, output_dir, scheme_dir):
    scheme = {}
    for line in open(files.joinpath(scheme_dir, "scheme.txt"), "r").read().splitlines():
        token = line.split("\t")
        scheme[token[0]] = token[1:]

    cols = []
    for filename in os.listdir(assemble_dir):
        assembly = filename.split(".")[0]
        source = files.joinpath(output_dir, "allele.{}.profile".format(assembly))
        d = (functional.seq(open(source, "r").read().splitlines())
             .map(lambda line: tuple(line.split("\t")))
             .to_dict())
        cols.append(pd.Series(d, name=assembly))
    table = pd.concat(cols, axis=1)

    for sch, loci in scheme.items():
        sink = files.joinpath(output_dir, "wgMLST_" + sch + ".tsv")
        table.loc[loci, :].to_csv(sink, sep="\t")


def profiling(output_dir, input_dir, db_dir, logger=None, threads=2, aligcov_cut=0.5, identity=90):
    if not logger:
        logger = logs.console_logger(__name__)

    logger.info("Renaming contigs...")
    assemble_dir = files.joinpath(output_dir, "query_assembly")
    files.create_if_not_exist(assemble_dir)
    namemap = rename(assemble_dir, input_dir)
    with open(files.joinpath(output_dir, "namemap.json"), "w") as f:
        f.write(json.dumps(namemap))

    logger.info("Profiling loci...")
    scheme_file = files.joinpath(db_dir, "scheme.csv")
    refseq_fna = files.joinpath(db_dir, "panRefSeq.fa")
    refseqdb_dir = files.joinpath(db_dir, "panRefSeq")
    blast_cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                  "sstart", "send", "evalue", "bitscore"]
    id_locus = lambda x, y: identify_locus(x, y, aligcov_cut, identity, blast_cols)
    profile_loci(refseq_fna, assemble_dir, output_dir, refseqdb_dir, id_locus, blast_cols)

    logger.info("Allocating alleles...")
    # blast_cols_str = "'6 {}'".format(" ".join(blast_cols))
    # bn = lambda x, y, z: NcbiblastnCommandline(query=x, db=y, out=z, num_threads=threads, outfmt=blast_cols_str)
    # allocate_alleles(assemble_dir, db_dir, output_dir, bn, blast_cols)

    logger.info("Collecting profiles by scheme...")
    # make_scheme_profile(assemble_dir, output_dir, scheme_dir)

    logger.info("Output...")
    # wgmlst_dir = files.joinpath(scheme_dir, "wgMLST")
    # files.clear_folder(wgmlst_dir)
    # os.system("mv {}/wgMLST_* {}".format(output_dir, wgmlst_dir))
    # os.system("rm -rf Assembly*")

