import os
import sys
import logging
import argparse
from pathlib import Path
from multiprocessing import Pool
from collections import Counter

import pandas as pd
from matplotlib import pyplot as plt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from BCBio import GFF

from utils import syscall, SQLiteDatabase, generate_record
from profiling import profiling, filter_alignments, sequence_alignment, sequence_encoder
plt.style.use('ggplot')


VERSION = {
    'roary': 'roary -w',
    'prokka': 'prokka -v 2>&1',
    'prodigal': 'prodigal -v 2>&1',
    'blastp': 'blastp -version | head -n 1',
}


class LoggerFactory:

    FMT = "%(asctime)-20s[%(levelname)s] %(message)s"
    DATEFMT = "%Y-%m-%d %H:%M:%S"

    def __init__(self):
        self._logger = logging.getLogger(__name__)
        self._logger.setLevel(logging.INFO)

    def addLogBoxHandler(self, logbox):
        log_handler = logging.Handler(logbox)
        formatter = logging.Formatter(self.FMT, self.DATEFMT)
        log_handler.setFormatter(formatter)
        self._logger.addHandler(log_handler)

    def addFileHandler(self, logfile):
        file_handler = logging.FileHandler(logfile, mode="a", encoding=None, delay=False)
        formatter = logging.Formatter(self.FMT, self.DATEFMT)
        file_handler.setFormatter(formatter)
        self._logger.addHandler(file_handler)

    def addConsoleHandler(self):
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(logging.INFO)
        formatter = logging.Formatter(self.FMT, self.DATEFMT)
        stream_handler.setFormatter(formatter)
        self._logger.addHandler(stream_handler)

    def create(self):
        return self._logger


def find_all_assemblies(assemblies):
    suffix = ('.fa', '.fna', '.fasta')
    return [x for x in assemblies if x.is_file() and x.suffix in suffix]


def clean_feature(src, dst):
    with open(src) as in_handle, open(dst, 'w') as out_handle:
        records = []
        for record in GFF.parse(in_handle):
            record.annotations = {}
            record.features = [feature for feature in record.features if feature.type == 'CDS']
            records.append(record)
        GFF.write(records, out_handle, include_fasta=True)


def annotate(assembly, output_dir, prodigaltf, locus_tag):
    if os.path.exists(output_dir):
        print(f'"Warning: Output {output_dir} already exists, may use old results."', file=sys.stderr)
        return
    cmd = ['prokka', '--prefix', 'prokka', '--cpus', '1', '--locustag', locus_tag, '--outdir', output_dir, assembly]
    if prodigaltf:
        cmd += ["--prodigaltf", prodigaltf]
    syscall(cmd, stdout=False)


def batch_process(function, args_ls, processes):
    p = Pool(processes=processes)
    for args in args_ls:
        p.apply_async(func=function, args=args)
    p.close()
    p.join()


def collect_gene_sequences(srcpath, outpath, gene_tag2id: dict):
    outpath.mkdir()
    for i in srcpath.iterdir():
        record_dict = SeqIO.to_dict(SeqIO.parse(i / 'prokka.ffn', 'fasta'))
        for gene_tag, record in record_dict.items():
            if gene_tag in gene_tag2id:
                with open(outpath / (gene_tag2id[gene_tag] + '.fa'), 'a') as handle:
                    handle.write(record.format('fasta'))


def plot_genome_coverage(data, figure_path, core=95):
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(data, bins=range(102), histtype='step', lw=2)
    ax.set_title("Genome coverage distribution", fontsize=22)
    ax.set_xlabel("Percentage of genomes covered by loci (%)", fontsize=18)
    ax.set_ylabel("Number of locus", fontsize=18)
    ax.set_xticks(list(range(0, 101, 10)) + [core])
    fig.savefig(figure_path, dpi=300, facecolor='w', bbox_inches='tight')


def parse_clustered_proteins(f):
    clusters = dict()
    with open(f) as handle:
        for line in handle:
            protein, cluster = line.split(': ')
            clusters[protein] = cluster.split()
    return clusters


def check_dependency():
    logger = logging.getLogger(__name__)
    for program_name, cmd in VERSION.items():
        child_process = syscall(cmd, stdout=True)
        if child_process.returncode:
            logger.error(msg=f"Could not determine version of {program_name}")
            sys.exit(0)
        else:
            # full_path = shutil.which(program_name)
            version = child_process.stdout.strip()
            logger.info(msg=f"Using {program_name:10} | {version}")


def main(input_files, output_path, threads, prodigaltf):
    logfile = os.path.join(output_path, 'benga.log')
    logger_factory = LoggerFactory()
    logger_factory.addFileHandler(logfile)
    logger_factory.addConsoleHandler()
    logger = logger_factory.create()
    check_dependency()
    all_assemblies = [Path(f) for f in input_files]
    output_path = Path(output_path)
    annot_dirname = output_path/'Annotated'
    annot_dirname.mkdir(exist_ok=True)

    all_assemblies = find_all_assemblies(all_assemblies)
    num_assemblies = len(all_assemblies)
    logger.info(f"Find {len(all_assemblies)} files")
    logger.info("Genomes annotation")
    args_ls = [(str(assembly), str(annot_dirname/assembly.stem), prodigaltf, assembly.stem) for assembly in all_assemblies]
    batch_process(annotate, args_ls, threads)

    gff_dirname = output_path/'GFF'
    gff_dirname.mkdir(exist_ok=True)
    for i in annot_dirname.iterdir():
        src = os.path.join(i, 'prokka.gff')
        dst = os.path.join(gff_dirname, i.name + '.gff')
        clean_feature(src, dst)

    logger.info("Gene cluster with roary")
    roary_path = output_path/'roary'
    syscall(f"roary -p {threads} -i 95 -s -f {roary_path} {gff_dirname/'*.gff'}")

    gene_presence_absence = pd.read_csv(roary_path/"gene_presence_absence.csv", index_col=0, usecols=range(13))
    gene_presence_absence.index = gene_presence_absence.index.str.replace(' ', '_')
    gene_frequency = dict(
        zip(gene_presence_absence.index, gene_presence_absence['No. isolates'].div(num_assemblies).mul(100).round(2))
    )
    gene_clusters = parse_clustered_proteins(os.path.join(roary_path, "clustered_proteins"))
    gene_tag2id = {gene_tag: gene_id for gene_id, gene_cluster in gene_clusters.items() for gene_tag in gene_cluster}

    logger.info(f"Collect gene sequences.")
    gene_sequences = roary_path / 'gene_sequences'
    collect_gene_sequences(annot_dirname, gene_sequences, gene_tag2id)
    rep_sequences = {}
    for i in gene_sequences.iterdir():
        rep_sequences[i.stem.replace(' ', '_')] = Counter((record.seq.ungap('-') for record in SeqIO.parse(i, 'fasta'))).most_common()[0][0]
    rep_records = [
        SeqRecord(gene_seq.translate(table=11, cds=True), id=gene_id)
        for gene_id, gene_seq in rep_sequences.items()
    ]
    filter_result = filter_alignments(sequence_alignment(rep_records, rep_records, threads))
    duplicated_genes = set()
    for qseqid, sseqid in filter_result:
        if gene_frequency[qseqid] > gene_frequency[sseqid]:
            duplicated_genes.add(sseqid)
        elif gene_frequency[qseqid] == gene_frequency[sseqid]:
            if qseqid.startswith('group'):
                duplicated_genes.add(qseqid)
            else:
                duplicated_genes.add(qseqid)
        else:
            duplicated_genes.add(sseqid)
    logger.info(f"Removo {len(duplicated_genes)} fragmented loci.")
    pangenome = gene_presence_absence.reindex(
        gene_presence_absence.index[~gene_presence_absence.index.isin(duplicated_genes)])
    pangenome['occurrence'] = pangenome['No. isolates'].div(num_assemblies).mul(100).round(2)
    pangenome['dna_seq'] = [str(rep_sequences[i]) for i in pangenome.index]

    pangenome['locus_id'] = [f'locus_{idx:05}' for idx, _ in enumerate(pangenome.index, 1)]

    database = os.path.join(output_path, 'sqlite3.db')
    sqlite_db = SQLiteDatabase(database)
    sqlite_db.create_table()
    sqlite_db.insert('scheme', zip(pangenome['locus_id'], pangenome['dna_seq']))
    sqlite_db.insert(
        'alleles',
        ((sequence_encoder(dna_seq), locus_id) for locus_id, dna_seq in zip(pangenome['locus_id'], pangenome['dna_seq']))
    )

    profile_path = output_path/'Profile'
    profile_path.mkdir(exist_ok=True)

    args_ls = []
    for i in annot_dirname.iterdir():
        seqfile = i/'prokka.fna'
        outfile = profile_path/(i.name + '.tsv')
        args_ls.append((seqfile, database, outfile, prodigaltf, 1))
    batch_process(profiling, args_ls, threads)

    counter = Counter()
    for x in profile_path.iterdir():
        profile = pd.read_csv(x, sep='\t')
        counter.update(profile.dropna()['locus_id'])

    pangenome['frequency'] = pangenome['locus_id'].map({locus_id: n for locus_id, n in counter.items()})
    pangenome['occurrence'] = pangenome['frequency'].div(num_assemblies).mul(100).round(2)
    logger.info(f"Removo {pangenome['frequency'].isna().sum()} inflate loci.")
    pangenome = pangenome[pangenome['frequency'].notna()]
    pangenome = pangenome[['Annotation', 'frequency', 'occurrence', 'dna_seq']]
    pangenome.to_csv(output_path/'pan_genome_info.txt', sep='\t', index=True)

    plot_genome_coverage(
        pangenome['occurrence'],
        output_path/'genome_coverage.png',
    )
    plot_genome_coverage(
        pangenome[pangenome['occurrence'] > 5]['occurrence'],
        output_path/'genome_coverage_5_prec.png',
    )
    logger.info(f"Down.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create wgMLST scheme")
    parser.add_argument("-i", "--input-files", nargs='+',
                        required=True,
                        help="Path of query assembly files.")
    parser.add_argument("-o", "--output-dir",
                        required=True,
                        help="Path of output directory.")
    parser.add_argument("--prodigaltf",
                        default='',
                        help="Path of prodigal training file. default:''")
    parser.add_argument("-t", "--threads",
                        type=int,
                        default=2,
                        help="Number of threads. default: 2")
    args = parser.parse_args()
    main(args.input_files, args.output_dir, args.threads, args.prodigaltf)
