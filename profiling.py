import os
import hashlib
import argparse
from tempfile import TemporaryDirectory
import pandas as pd
from Bio import SeqIO
from utils import syscall, make_blast_database, run_blastp, SQLiteDatabase, generate_record


def sequence_encoder(sequence):
    """
    Convert nucleotide sequence to hash codes by hash function 'sha256'.
    """
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def find_genes(input_fasta, training_file=None):
    """
    Protein-coding gene prediction.
    """
    records = []
    with TemporaryDirectory() as tmpdir:
        prodigal_output = os.path.join(tmpdir, 'genes.fna')
        cmd = ['prodigal', '-i', input_fasta, '-d', prodigal_output, '-c', '-m', '-q']
        if training_file:
            cmd += ['-t', training_file]
        syscall(cmd)
        for record in SeqIO.parse(prodigal_output, 'fasta'):
            record.id = sequence_encoder(str(record.seq))
            record.seq = record.seq.translate(table=11)
            records.append(record)
    return records


def sequence_alignment(query_records, subject_records, threads):
    with TemporaryDirectory() as tmpdir:
        query, subject = os.path.join(tmpdir, 'query.faa'), os.path.join(tmpdir, 'subject.faa')
        SeqIO.write(query_records, query, 'fasta')
        SeqIO.write(subject_records, subject, 'fasta')
        make_blast_database(subject, subject)
        stdout, stderr = run_blastp(query, subject, threads)
    return stdout.splitlines()


def filter_alignments(aligns, identity=95, min_cov=.75, max_cov=1.25):
    """
    Filter identity >= 95, ratio of query sequence length to subject sequence length between 0.75 and 1.25
    """
    result = []
    for line in aligns:
        qseqid, sseqid, pident, length, qlen, slen = line.strip().split()
        pident, length, qlen, slen = float(pident), float(length), float(qlen), float(slen)
        if qseqid != sseqid and pident >= identity and min_cov <= qlen/slen < max_cov and min_cov <= qlen/length < max_cov:
            result.append((qseqid, sseqid))
    return result


def profiling(infile, database, outfile, training_file, threads, identity=95):
    sqlite_db = SQLiteDatabase(database)

    gene_records = find_genes(infile, training_file)
    existed_alleles = sqlite_db.search([gene.id for gene in gene_records])
    existed_allele_ids = set(x[0] for x in existed_alleles)
    candidates_records = [record for record in gene_records if record.id not in existed_allele_ids]
    scheme_records = [generate_record(locus_tag, dna_seq, True) for locus_tag, dna_seq in sqlite_db.fetch_scheme()]
    new_alleles = filter_alignments(
        sequence_alignment(candidates_records, scheme_records, threads), identity=identity
    )

    if new_alleles:
        sqlite_db.insert('alleles', new_alleles)
    df = pd.DataFrame(existed_alleles + new_alleles, columns=['allele_id', 'locus_id'])
    df = df.sort_values('allele_id', kind='mergesort').drop_duplicates('locus_id')
    df = df.set_index('locus_id').reindex((x.id for x in scheme_records)).sort_index()
    df.to_csv(outfile, sep='\t')


def main():
    parser = argparse.ArgumentParser('Benga')
    parser.add_argument("-i", "--input",
                        required=True,
                        help="Path of query genome.")
    parser.add_argument("-o", "--output",
                        required=True,
                        help="Path of output file.")
    parser.add_argument("-d", "--database",
                        required=True,
                        help="Path of core-genome MLST database.")
    parser.add_argument("--prodigaltf",
                        default='',
                        help="Path of prodigal training file. default: None")
    parser.add_argument("-t", "--threads",
                        type=int,
                        default=1,
                        help="Number of threads. default: 1")
    args = parser.parse_args()

    profiling(
        infile=args.input,
        database=args.database,
        outfile=args.output,
        training_file=args.prodigaltf,
        threads=args.threads,
    )


if __name__ == '__main__':
    main()
