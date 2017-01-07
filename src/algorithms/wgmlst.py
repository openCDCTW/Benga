import json
import os
import shutil
from collections import defaultdict, OrderedDict

import functional
import pandas as pd
from Bio import SeqIO

from src.utils import files, seq


def renumber(assemble_dir, query_dir):
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


def profile_loci(blastquery, assemble_dir, output_dir, qb):
    seqlen = (functional.seq(SeqIO.parse(blastquery, "fasta"))
              .map(lambda rec: (rec.id, len(rec.seq)))
              .to_dict())
    for filename in os.listdir(assemble_dir):  # TODO: parallelize
        m = (functional.seq(seqlen.keys())
             .map(lambda id: (id, ""))
             .sorted(key=lambda x: x[0])
             .to_list())
        blast = qb(OrderedDict(m), output_dir, files.joinpath(assemble_dir, filename), blastquery, seqlen)
        output = (functional.seq(blast.items())
                  .map(lambda x: (x[0], 1 if len(x[1]) > 1 else 0))
                  .map(lambda x: x[0] + "\t" + str(x[1]))
                  .to_list())

        name = filename.split(".")[0]
        with open(files.joinpath(output_dir, "locusAP." + name + ".list"), "w") as file:
            file.write("\n".join(output))


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
        os.system(makeblastdb(blastnfile, "nucl", files.joinpath(output_dir, "AssemblyDB_" + qfn)))

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


def makeblastdb(input, dbtype, output):
    return "makeblastdb -in {} -dbtype {} -out {}".format(input, dbtype, output)


def blastn(num_threads, query, db, output, outfmt):
    return "blastn -num_threads {} -query {} -db {} -out {} -outfmt {}".format(num_threads, query, db, output, outfmt)


def query_blast(m, dest_dir, dbsource, query, seqlen, aligcov_cut, pident_cut, threads, cols):
    blast_out = files.joinpath(dest_dir, "blast.out")
    locidb = files.joinpath(dest_dir, "lociDB")

    os.system(makeblastdb(dbsource, "nucl", locidb))
    os.system(blastn(threads, query, locidb, blast_out, "'6 {items}'".format(items=" ".join(cols))))

    result = pd.read_csv(blast_out, sep="\t", header=None, names=cols)
    result["qlen"] = list(map(lambda x: seqlen[x], result["qseqid"]))
    result["aligcov"] = (result["length"] - result["gapopen"]) / result["qlen"]
    result = result[(result["aligcov"] >= aligcov_cut) & (result["pident"] >= pident_cut)]
    for index, row in result.iterrows():
        row_elements = [row["aligcov"], row["length"], row["qlen"], row["pident"], row["sseqid"], row["sstart"],
                        row["send"]]
        m[row["qseqid"]] = "\t".join(map(str, row_elements))
    os.remove(blast_out)
    return m


def profiling(output_dir, input_dir, db_dir, logger=None, blastn_threads=2, aligcov_cut=0.5, pident_cut=90):
    if logger:
        logger.info("Reformating contigs...")
    assemble_dir = files.joinpath(output_dir, "query_assembly")
    files.create_if_not_exist(assemble_dir)
    namemap = renumber(assemble_dir, input_dir)
    with open(files.joinpath(output_dir, "namemap.json"), "w") as f:
        f.write(json.dumps(namemap))

    if logger:
        logger.info("Profiling loci...")
    scheme_dir = files.joinpath(output_dir, "scheme")
    shutil.copytree(files.joinpath(db_dir, "scheme"), scheme_dir)
    blastquery = files.joinpath(scheme_dir, "pan.txt")
    blast_cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                  "sstart", "send", "evalue", "bitscore"]
    qb = lambda a, b, c, d, e: query_blast(a, b, c, d, e, aligcov_cut, pident_cut, blastn_threads, blast_cols)
    profile_loci(blastquery, assemble_dir, output_dir, qb)

    if logger:
        logger.info("Allocating alleles...")
    blast_cols_str = "'6 {}'".format(" ".join(blast_cols))
    bn = lambda x, y, z: blastn(blastn_threads, x, y, z, blast_cols_str)
    allocate_alleles(assemble_dir, db_dir, output_dir, bn, blast_cols)

    if logger:
        logger.info("Collecting profiles by scheme...")
    make_scheme_profile(assemble_dir, output_dir, scheme_dir)

    if logger:
        logger.info("Output...")
    wgmlst_dir = files.joinpath(scheme_dir, "wgMLST")
    files.clear_folder(wgmlst_dir)
    os.system("mv {}/wgMLST_* {}".format(output_dir, wgmlst_dir))
    os.system("rm -rf Assembly*")