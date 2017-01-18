import argparse
import json
import os
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
    mapping = matrix[matrix["No. isolates"] == matrix["No. sequences"]]
    mapping.to_csv(files.joinpath(output_dir, "locusmapping.txt"), sep="\t")

    paralog = matrix[matrix["No. isolates"] != matrix["No. sequences"]]
    paralog.to_csv(files.joinpath(output_dir, "paralog.txt"), sep="\t", index=False)


def parse_csv(output_dir, ffn_dir, locus_dir):
    roary_result = pd.read_csv(files.joinpath(output_dir, "roary", "gene_presence_absence.csv"))
    roary_result.index = pd.Index(map(lambda x: "SAL{0:07d}".format(x + 1), roary_result.index), name="locus")

    matrix = roary_result[["Gene", "No. isolates", "No. sequences", "Annotation"]]
    save_locus_map(matrix, output_dir)

    for colname, col in roary_result.iloc[:, 14:].iteritems():
        handle = SeqIO.parse(files.joinpath(ffn_dir, colname + ".ffn"), "fasta")
        cds = {record.id: record.seq for record in handle}
        for rowname, loci in col.iteritems():
            if isinstance(loci, str):
                with open(files.joinpath(locus_dir, rowname + ".tmp"), "a") as file:
                    # write single-line FASTA for fastx_collapser, which doesn't accept multi-line
                    prokkaID = loci.split("___")[0]
                    file.write(">" + prokkaID + "\n" + str(cds[prokkaID]) + "\n")


def find_longest_seq(filename):
    return (functional.seq(SeqIO.parse(filename, "fasta"))
            .map(lambda record: str(record.seq))
            .max_by(lambda s: len(s)))


def pan_refseq(database_dir, locus_files, locus_dir):
    """find longest sequence in a locus file as the pan RefSeq"""
    refseqs = (functional.seq(locus_files)
               .map(lambda file: os.path.splitext(file)[0])
               .map(lambda locus: (locus, files.joinpath(locus_dir, locus + ".fa")))
               .map(lambda file: (file[0], find_longest_seq(file[1])))
               .to_dict())

    records = [seq.new_record(key, value) for key, value in refseqs.items()]
    SeqIO.write(records, files.joinpath(database_dir, "panRefSeq.fa"), "fasta")

    with open(files.joinpath(database_dir, "panRefSeq.json"), "w") as file:
        file.write(json.dumps(refseqs))
    return refseqs


def dispatch_loci(mapping_file, total_isolates):
    dispatcher = defaultdict(list)
    mapping = pd.read_csv(mapping_file, sep="\t")
    for index, row in mapping.iterrows():
        locus = row["locus"]
        occ = row["No. isolates"]
        occrate = occ / total_isolates

        dispatcher["pan"].append(locus)
        if 0.95 <= occrate <= 1:
            dispatcher["core"].append(locus)
    return dispatcher


def save_schemes(dispatcher, refseqs, scheme_dir):
    schemes = []
    for scheme, loci in dispatcher.items():
        if len(loci) != 0:
            records = [seq.new_record(l, refseqs[l]) for l in loci]
            SeqIO.write(records, files.joinpath(scheme_dir, scheme + ".txt"), "fasta")
            content = "\t".join(loci)
            schemes.append(scheme + "\t" + content)
    with open(files.joinpath(scheme_dir, "scheme.txt"), "w") as file:
        file.write("\n".join(schemes))


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
        cmd.roary(output_dir, output_dir, threads, ident_min)

    if logger:
        logger.info("Parsing csv...")
    locus_dir = files.joinpath(output_dir, "locusfiles")
    files.create_if_not_exist(locus_dir)
    parse_csv(output_dir, ffn_dir, locus_dir)

    if logger:
        logger.info("Reducing identical sequences in a FASTA file...")
    locus_files = [file for file in os.listdir(locus_dir) if file.endswith(".tmp")]
    if use_docker:
        docker.fastx(locus_files, output_dir)  # TODO: bug -- docker never stops
    else:
        c = (functional.seq(locus_files)
             .map(lambda f: (f, os.path.splitext(f)[0]))
             .map(lambda x: cmd.fastx(files.joinpath(locus_dir, x[0]),
                                      files.joinpath(locus_dir, x[1] + ".fa")))
             .to_list())
        with ProcessPoolExecutor() as executor:
            executor.map(os.system, c)
        for file in locus_files:
            os.remove(files.joinpath(locus_dir, file))

    if logger:
        logger.info("Find longest sequence as pan RefSeq...")
    refseqs = pan_refseq(database_dir, locus_files, locus_dir)

    if logger:
        logger.info("Making schemes...")
    scheme_dir = files.joinpath(output_dir, "scheme")
    files.create_if_not_exist(scheme_dir)

    total_isolates = len(os.listdir(ffn_dir))
    mapping_file = files.joinpath(output_dir, "locusmapping.txt")
    dispatcher = dispatch_loci(mapping_file, total_isolates)
    save_schemes(dispatcher, refseqs, scheme_dir)

    if logger:
        logger.info("Collecting outputs...")
    shutil.copy(files.joinpath(output_dir, "roary", "summary_statistics.txt"), database_dir)
    shutil.copy(files.joinpath(output_dir, "locusmapping.txt"), database_dir)
    shutil.move(scheme_dir, database_dir)
    for file in locus_files:
        if file.endswith(".tmp"):
            os.remove(files.joinpath(locus_dir, file))
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
