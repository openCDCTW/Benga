import subprocess

from Bio.Alphabet import generic_dna, generic_protein
from Bio.Blast.Applications import NcbiblastnCommandline, NcbiblastpCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import write

BLAST_COLUMNS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                 "sstart", "send", "evalue", "bitscore"]


def new_record(seqid, seq, desc="", seqtype=None):
    seqtype = generic_protein if seqtype == "protein" else generic_dna
    if type(seq) == str:
        return SeqRecord(Seq(seq, seqtype), id=seqid, description=desc)
    elif type(seq) == Seq:
        return SeqRecord(seq, id=seqid, description=desc)
    else:
        print("None supported type: {}".format(type(seq)))


def replace_id(record, newid):
    return new_record(newid, str(record.seq))


def save_records(seqs, filename):
    write(seqs, filename, "fasta")


def compile_blastndb(input_file, output_file):
    subprocess.run(["makeblastdb", "-in", input_file, "-dbtype", "nucl", "-out", output_file], shell=True)


def compile_blastpdb(input_file, output_file):
    subprocess.run(["makeblastdb", "-in", input_file, "-dbtype", "prot", "-out", output_file], shell=True)


def query_blastndb(query, db_dir, output_file, cols, threads=2):
    NcbiblastnCommandline(query=query, db=db_dir, out=output_file,
                          outfmt="'6 {}'".format(" ".join(cols)), num_threads=threads)()


def query_blastpdb(query, db_dir, output_file, cols, cov=90, threads=2):
    NcbiblastpCommandline(query=query, db=db_dir, out=output_file,
                          outfmt="'6 {}'".format(" ".join(cols)), num_threads=threads,
                          qcov_hsp_perc=cov)()
