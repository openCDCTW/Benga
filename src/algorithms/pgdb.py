import argparse
import json
import os
import heapq
import shutil
from collections import defaultdict
import functional
import pandas as pd
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor

from src.utils import files, seq, docker, cmd


def parse_filenames(path, ext=".fna"):
    return [name for name in os.listdir(path) if name.endswith(ext)]


def format_contgis(filenames, input_dir, working_dir):
    namemap = {}
    for i, oldname in enumerate(filenames, 1):
        newname = "Assembly_" + str(i) + ".fa"
        namemap[oldname] = newname
        with open(files.joinpath(working_dir, newname), "w") as file:
            for j, contig in enumerate(SeqIO.parse(files.joinpath(input_dir, oldname), "fasta"), 1):
                seqid = "A_{i}::C_{j}".format(**locals())
                SeqIO.write(seq.replace_id(contig, seqid), file, "fasta")
    return namemap


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


def save_locus_map(matrix, output_dir):
    mapping_file = files.joinpath(output_dir, "locusmapping.txt")
    matrix[matrix["No. isolates"] == matrix["No. sequences"]].to_csv(mapping_file, sep="\t")

    paralog_file = files.joinpath(output_dir, "paralog.txt")
    matrix[matrix["No. isolates"] != matrix["No. sequences"]].to_csv(paralog_file, sep="\t")


def dump_alleles(row, filename):
    records = [seq.new_record(str(i), a) for i, a in enumerate(row, 1)]
    seq.save_records(records, filename)


def longest_alleles(row):
    alleles = list(map(lambda x: (len(x), x), row))
    return heapq.nlargest(1, alleles)[0][1]


def parse_csv(output_dir, ffn_dir, locus_dir, metadata_colnumber=14):
    roary_result = pd.read_csv(files.joinpath(output_dir, "roary", "gene_presence_absence.csv"))
    roary_result.index = pd.Index(map(lambda x: "SAL{0:07d}".format(x + 1), roary_result.index), name="locus")
    total_isolates = len(roary_result.columns) - metadata_colnumber

    matrix = roary_result[["Gene", "No. isolates", "No. sequences", "Annotation"]]
    save_locus_map(matrix, output_dir)

    matrix = roary_result[roary_result["No. isolates"] == roary_result["No. sequences"]]
    allele_map = defaultdict(set)
    for colname, col in matrix.iloc[:, metadata_colnumber:].iteritems():
        ffn_file = files.joinpath(ffn_dir, colname + ".ffn")
        cds = {record.id: record.seq for record in SeqIO.parse(ffn_file, "fasta")}
        for locus, allele in col.iteritems():
            if type(allele) == str:
                allele_map[locus].add(str(cds[allele]))

    longests = []
    for locus, alleles in allele_map.items():
        dump_alleles(alleles, files.joinpath(locus_dir, locus + ".fa"))
        longests[locus] = longest_alleles(alleles)
    return longests, total_isolates


def make_schemes(mapping_file, refseqs, total_isolates, scheme_dir):
    mapping = pd.read_csv(mapping_file, sep="\t")
    mapping["occurence"] = list(map(lambda x: round(x/total_isolates * 100, 2), mapping["No. isolates"]))
    mapping["sequence"] = list(map(lambda x: refseqs[x], mapping["locus"]))
    mapping[["locus", "occurence", "sequence"]].to_csv(files.joinpath(scheme_dir, "scheme.csv"), index=False)


def create_new_locusfiles(database_dir):
    locusfiles_dir = files.joinpath(database_dir, "locusfiles")
    for filename in os.listdir(locusfiles_dir):
        filepath = files.joinpath(locusfiles_dir, filename)
        records = [seq.new_record(filename + "::" + rec.id, str(rec.seq)) for rec in SeqIO.parse(filepath, "fasta")]
        SeqIO.write(records, files.joinpath(locusfiles_dir, filename + ".new"), "fasta")


def annotate_configs(input_dir, output_dir, logger=None, use_docker=True):
    if logger:
        logger.info("Formating contigs...")
    filenames = parse_filenames(input_dir)

    assembly_dir = files.joinpath(output_dir, "Assembly")
    files.create_if_not_exist(assembly_dir)
    namemap = format_contgis(filenames, input_dir, assembly_dir)

    if logger:
        logger.info("Annotating...")
    annotate_dir = files.joinpath(output_dir, "Annotated")
    files.create_if_not_exist(annotate_dir)
    if use_docker:
        docker.prokka(namemap.values(), annotate_dir, assembly_dir)
    else:
        c = [cmd.prokka(x, annotate_dir, assembly_dir) for x in namemap.values()]
        with ProcessPoolExecutor() as executor:
            executor.map(os.system, c)

    # protein CDS file
    if logger:
        logger.info("Moving ffn files...")
    ffn_dir = files.joinpath(output_dir, "FFN")
    files.create_if_not_exist(ffn_dir)
    move_file(annotate_dir, ffn_dir, ".ffn")

    # annotation file
    if logger:
        logger.info("Moving gff files...")
    gff_dir = files.joinpath(output_dir, "GFF")
    files.create_if_not_exist(gff_dir)
    move_file(annotate_dir, gff_dir, ".gff")

    if logger:
        logger.info("Creating nonCDS.json...")
    create_noncds(output_dir, gff_dir)


def make_database(output_dir, logger=None, threads=2, ident_min=95, use_docker=True):
    database_dir = files.joinpath(output_dir, "DB")
    files.create_if_not_exist(database_dir)
    ffn_dir = files.joinpath(output_dir, "FFN")

    if logger:
        logger.info("Calculating the pan genome")
    if use_docker:
        docker.roary(output_dir, threads=threads, ident_min=ident_min)
    else:
        c = cmd.roary(output_dir, files.joinpath(output_dir, "GFF"), threads, ident_min)
        os.system(c)

    if logger:
        logger.info("Parsing csv and finding longest allele as RefSeq...")
    locus_dir = files.joinpath(output_dir, "locusfiles")
    files.create_if_not_exist(locus_dir)
    refseqs, total_isolates = parse_csv(output_dir, ffn_dir, locus_dir)

    if logger:
        logger.info("Saving pan RefSeq...")
    records = [seq.new_record(key, str(value)) for key, value in refseqs.items()]
    SeqIO.write(records, files.joinpath(database_dir, "panRefSeq.fa"), "fasta")

    if logger:
        logger.info("Making dynamic schemes...")
    mapping_file = files.joinpath(output_dir, "locusmapping.txt")
    make_schemes(mapping_file, refseqs, total_isolates, database_dir)

    if logger:
        logger.info("Collecting outputs...")
    shutil.copy(files.joinpath(output_dir, "roary", "summary_statistics.txt"), database_dir)
    shutil.copy(files.joinpath(output_dir, "locusmapping.txt"), database_dir)
    shutil.move(locus_dir, database_dir)

    create_new_locusfiles(database_dir)

    if logger:
        logger.info("Done!!")


def parse_args():
    arg_parser = argparse.ArgumentParser()

    arg_parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output data directory. (necessary)"
    )

    arg_parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input data directory. (necessary)"
    )

    arg_parser.add_argument(
        "-t", "--threads",
        type=int,
        default=2,
        help="Number of threads for computation. [Default: 2]",
        metavar="FULL_FILE_PATH"
    )

    arg_parser.add_argument(
        "--identity",
        help="The minimum percentage identity for blastp. [Default: 95]",
        type=int,
        default=95,
        metavar='THRESHOLD'
    )

    return arg_parser.parse_args()


def main():
    args = parse_args()

    input_dir = args.input
    output_dir = args.output
    threads = args.threads
    ident_min = args.identity

    annotate_configs(input_dir, output_dir)
    make_database(output_dir, threads=threads, ident_min=ident_min)


if __name__ == "__main__":
    main()
