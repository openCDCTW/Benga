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
        newname = "Genome_{}.fa".format(i)
        namemap[oldname] = newname
        with open(files.joinpath(working_dir, newname), "w") as file:
            for j, contig in enumerate(SeqIO.parse(files.joinpath(input_dir, oldname), "fasta"), 1):
                seqid = "G_{}::C_{}".format(i, j)
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


def extract_profiles(roary_matrix_file, locusmeta_file, paralogmeta_file, metadata_cols=14):
    matrix = pd.read_csv(roary_matrix_file)
    matrix.index = pd.Index(map(lambda x: "SAL{0:07d}".format(x + 1), matrix.index), name="locus")

    appear_once_locus = matrix[matrix["No. isolates"] == matrix["No. sequences"]]
    not_appear_once_locus = matrix[matrix["No. isolates"] != matrix["No. sequences"]]
    save_metadata(appear_once_locus, locusmeta_file)
    save_metadata(not_appear_once_locus, paralogmeta_file)

    profiles = appear_once_locus.iloc[:, metadata_cols:]
    isolates = len(matrix.columns) - metadata_cols
    return profiles, isolates


def save_metadata(loci, meta_file, select_col=None):
    if not select_col:
        select_col = ["Gene", "No. isolates", "No. sequences", "Annotation"]
    loci[select_col].to_csv(meta_file, sep="\t")


def collect_record_profile(args):
    subject, profile, ffn_dir = args
    ffn_file = files.joinpath(ffn_dir, "{}.ffn".format(subject))
    seqs = {record.id: record.seq for record in SeqIO.parse(ffn_file, "fasta")}
    new_profile = pd.Series(name=subject)
    for locus, allele in profile.dropna().iteritems():
        s = str(seqs[allele])
        rec = seq.new_record(operations.make_seqid(s), s)
        new_profile.set_value(locus, rec)
    return new_profile


def generate_record_profiles(profiles, ffn_dir, threads):
    args = ((subject, profile, ffn_dir) for subject, profile in profiles.iteritems())
    with ProcessPoolExecutor(threads) as executor:
        new_profiles = list(executor.map(collect_record_profile, args))
    return pd.concat(new_profiles, axis=1).sort_index().sort_index(axis=1)


def allele_infomation(args):
    locus, alleles, locus_dir = args
    alleles = [(rec.id, str(rec.seq)) for rec in alleles.dropna().values]
    freq = Counter(alleles)
    refseq = freq.most_common(1)[0][0][1]
    frequency = {allele[0]: count for allele, count in freq.items()}
    records = [seq.new_record(allele[0], allele[1]) for allele in freq.keys()]
    seq.save_records(records, files.joinpath(locus_dir, locus + ".fa"))
    return locus, frequency, refseq


def generate_allele_infos(record_profiles, allele_freq_file, locus_dir, refseq_file, threads):
    frequency = {}
    refseqs = {}
    args = [(locus, alleles, locus_dir) for locus, alleles in record_profiles.iterrows()]
    with ProcessPoolExecutor(threads) as executor:
        for locus, freq, refseq in executor.map(allele_infomation, args):
            frequency[locus] = freq
            refseqs[locus] = refseq

    with open(allele_freq_file, "w") as file:
        file.write(json.dumps(frequency))

    records = [seq.new_record(key, str(value)) for key, value in refseqs.items()]
    SeqIO.write(records, refseq_file, "fasta")
    return refseqs


def make_schemes(locusmeta_file, scheme_file, refseqs, total_isolates):
    mapping = pd.read_csv(locusmeta_file, sep="\t")
    mapping["occurence"] = list(map(lambda x: round(x/total_isolates * 100, 2), mapping["No. isolates"]))
    mapping["sequence"] = list(map(lambda x: refseqs[x], mapping["locus"]))
    mapping[["locus", "occurence", "sequence"]].to_csv(scheme_file, index=False, sep="\t")


def annotate_configs(input_dir, output_dir, logger=None, threads=8, use_docker=True):
    if not logger:
        logger = logs.console_logger(__name__)

    logger.info("Formating contigs...")
    filenames = parse_filenames(input_dir)

    genome_dir = files.joinpath(output_dir, "Genomes")
    files.create_if_not_exist(genome_dir)
    namemap = format_contigs(filenames, input_dir, genome_dir)
    with open(files.joinpath(output_dir, "namemap.json"), "w") as f:
        f.write(json.dumps(namemap))

    logger.info("Annotating...")
    annotate_dir = files.joinpath(output_dir, "Annotated")
    files.create_if_not_exist(annotate_dir)
    if use_docker:
        docker.prokka(genome_dir, annotate_dir)
    else:
        c = [cmds.form_prokka_cmd(x, genome_dir, annotate_dir) for x in namemap.values()]
        with ProcessPoolExecutor(int(threads / 2)) as executor:
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

    database_dir = files.joinpath(output_dir, "database")
    files.create_if_not_exist(database_dir)

    logger.info("Calculating the pan genome...")
    min_identity = 95
    if use_docker:
        docker.roary(files.joinpath(output_dir, "GFF"), output_dir, min_identity, threads)
    else:
        c = cmds.form_roary_cmd(files.joinpath(output_dir, "GFF"), output_dir, min_identity, threads)
        os.system(c)

    logger.info("Extract profiles from roary result matrix...")
    matrix_file = files.joinpath(output_dir, "roary", "gene_presence_absence.csv")
    locusmeta_file = files.joinpath(database_dir, "locus_metadata.tsv")
    paralogmeta_file = files.joinpath(database_dir, "paralog_metadata.tsv")
    profiles, total_isolates = extract_profiles(matrix_file, locusmeta_file, paralogmeta_file)

    logger.info("Generating allele profiles...")
    ffn_dir = files.joinpath(output_dir, "FFN")
    profile_file = files.joinpath(database_dir, "allele_profiles.tsv")
    record_profiles = generate_record_profiles(profiles, ffn_dir, threads)
    record_profiles.applymap(lambda x: x.id if type(x) != float else x).to_csv(profile_file, sep="\t")

    logger.info("Generating allele sequences, allele frequencies and reference sequence...")
    refseq_file = files.joinpath(database_dir, "panRefSeq.fa")
    allele_freq_file = files.joinpath(database_dir, "allele_frequency.json")
    locus_dir = files.joinpath(database_dir, "locusfiles")
    files.create_if_not_exist(locus_dir)
    refseqs = generate_allele_infos(record_profiles, allele_freq_file, locus_dir, refseq_file, threads)

    logger.info("Making dynamic schemes...")
    scheme_file = files.joinpath(database_dir, "scheme.tsv")
    make_schemes(locusmeta_file, scheme_file, refseqs, total_isolates)
    logger.info("Done!!")

