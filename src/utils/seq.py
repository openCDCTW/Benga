from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import write


def new_record(seqid, seq, desc=""):
    return SeqRecord(Seq(seq, generic_dna), id=seqid, description=desc)


def replace_id(record, newid):
    return new_record(newid, str(record.seq))


def save_records(seqs, filename):
    write(seqs, filename, "fasta")
