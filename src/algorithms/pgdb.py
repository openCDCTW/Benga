import json
import os
import shutil
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from Bio import SeqIO

from src.models import logs
from src.utils import files, seq, docker, cmds
from src.utils import operations


def parse_filenames(path, ext=".fna"):
    return [name for name in os.listdir(path) if name.endswith(ext)]


def format_contigs(filenames, input_dir, working_dir):
    namemap = {}
    for i, oldname in enumerate(filenames, 1):
        newname = "Assembly_{}.fa".format(i)
        namemap[oldname] = newname
        with open(files.joinpath(working_dir, newname), "w") as file:
            for j, contig in enumerate(SeqIO.parse(files.joinpath(input_dir, oldname), "fasta"), 1):
                seqid = "A_{}::C_{}".format(i, j)
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


def dump_alleles(row, filename):
    records = [seq.new_record(operations.make_seqid(x), x) for x in row]
    seq.save_records(records, filename)


def most_frequent_allele(row):
    counter = Counter(row)
    return counter.most_common(1)[0][0]


def identify_pan_refseq(output_dir, ffn_dir, locus_dir):
    profiles, isolates = extract_profiles(output_dir)
    new_profiles = pd.DataFrame(columns=profiles.columns, index=profiles.index)
    allele_map = defaultdict(set)
    for colname, col in profiles.iteritems():
        ffn_file = files.joinpath(ffn_dir, colname + ".ffn")
        cds = {record.id: record.seq for record in SeqIO.parse(ffn_file, "fasta")}
        for locus, allele in col.iteritems():
            if type(allele) == str:
                s = str(cds[allele])
                allele_map[locus].add(s)
                new_profiles.loc[locus, colname] = operations.make_seqid(s)
    new_profiles.to_csv(files.joinpath(output_dir, "allele_profiles.tsv"), sep="\t")

    frequent = {}
    for locus, alleles in allele_map.items():
        dump_alleles(alleles, files.joinpath(locus_dir, locus + ".fa"))
        frequent[locus] = most_frequent_allele(alleles)
    return frequent, isolates


def extract_profiles(output_dir, metadata_colnumber=14):
    profile_matrix = pd.read_csv(files.joinpath(output_dir, "roary", "gene_presence_absence.csv"))
    profile_matrix.index = pd.Index(map(lambda x: "SAL{0:07d}".format(x + 1), profile_matrix.index), name="locus")
    isolates = len(profile_matrix.columns) - metadata_colnumber
    appear_once_locus = divide_matrix(profile_matrix, output_dir)
    profiles = appear_once_locus.iloc[:, metadata_colnumber:]
    return profiles, isolates


def divide_matrix(mat, output_dir):
    appear_once_locus = mat[mat["No. isolates"] == mat["No. sequences"]]
    not_appear_once_locus = mat[mat["No. isolates"] != mat["No. sequences"]]
    save_metadata(appear_once_locus, not_appear_once_locus, output_dir)
    return appear_once_locus


def save_metadata(appear_once_locus, not_appear_once_locus, output_dir):
    select_col = ["Gene", "No. isolates", "No. sequences", "Annotation"]
    metadata = appear_once_locus[select_col]
    paralog = not_appear_once_locus[select_col]
    metadata.to_csv(files.joinpath(output_dir, "locus_metadata.tsv"), sep="\t")
    paralog.to_csv(files.joinpath(output_dir, "paralog_metadata.tsv"), sep="\t")


def make_schemes(mapping_file, refseqs, total_isolates, scheme_dir):
    mapping = pd.read_csv(mapping_file, sep="\t")
    mapping["occurence"] = list(map(lambda x: round(x/total_isolates * 100, 2), mapping["No. isolates"]))
    mapping["sequence"] = list(map(lambda x: refseqs[x], mapping["locus"]))
    mapping[["locus", "occurence", "sequence"]].to_csv(files.joinpath(scheme_dir, "scheme.csv"), index=False)


def annotate_configs(input_dir, output_dir, logger=None, use_docker=True):
    if not logger:
        logger = logs.console_logger(__name__)

    logger.info("Formating contigs...")
    filenames = parse_filenames(input_dir)

    assembly_dir = files.joinpath(output_dir, "Assembly")
    files.create_if_not_exist(assembly_dir)
    namemap = format_contigs(filenames, input_dir, assembly_dir)
    with open(files.joinpath(output_dir, "namemap.json"), "w") as f:
        f.write(json.dumps(namemap))

    logger.info("Annotating...")
    annotate_dir = files.joinpath(output_dir, "Annotated")
    files.create_if_not_exist(annotate_dir)
    if use_docker:
        docker.prokka(assembly_dir, annotate_dir)
    else:
        c = [cmds.form_prokka_cmd(x, assembly_dir, annotate_dir) for x in namemap.values()]
        with ProcessPoolExecutor() as executor:
            executor.map(os.system, c)

    logger.info("Moving protein CDS (.ffn) files...")
    ffn_dir = files.joinpath(output_dir, "FFN")
    files.create_if_not_exist(ffn_dir)
    move_file(annotate_dir, ffn_dir, ".ffn")

    logger.info("Moving annotation (.gff) files...")
    gff_dir = files.joinpath(output_dir, "GFF")
    files.create_if_not_exist(gff_dir)
    move_file(annotate_dir, gff_dir, ".gff")

    logger.info("Creating nonCDS.json...")
    create_noncds(output_dir, gff_dir)


def make_database(output_dir, logger=None, threads=2, use_docker=True):
    if not logger:
        logger = logs.console_logger(__name__)

    database_dir = files.joinpath(output_dir, "DB")
    files.create_if_not_exist(database_dir)

    logger.info("Calculating the pan genome...")
    min_identity = 95
    if use_docker:
        docker.roary(files.joinpath(output_dir, "GFF"), output_dir, min_identity, threads)
    else:
        c = cmds.form_roary_cmd(files.joinpath(output_dir, "GFF"), output_dir, min_identity, threads)
        os.system(c)

    logger.info("Finding most frequent allele as RefSeq...")
    ffn_dir = files.joinpath(output_dir, "FFN")
    locus_dir = files.joinpath(output_dir, "locusfiles")
    files.create_if_not_exist(locus_dir)
    refseqs, total_isolates = identify_pan_refseq(output_dir, ffn_dir, locus_dir)

    logger.info("Saving pan RefSeq...")
    records = [seq.new_record(key, str(value)) for key, value in refseqs.items()]
    SeqIO.write(records, files.joinpath(database_dir, "panRefSeq.fa"), "fasta")

    logger.info("Making dynamic schemes...")
    mapping_file = files.joinpath(output_dir, "locus_metadata.tsv")
    make_schemes(mapping_file, refseqs, total_isolates, database_dir)

    logger.info("Collecting outputs...")
    shutil.copy(files.joinpath(output_dir, "roary", "summary_statistics.txt"), database_dir)
    shutil.copy(files.joinpath(output_dir, "locus_metadata.tsv"), database_dir)
    shutil.move(locus_dir, database_dir)
    logger.info("Done!!")

