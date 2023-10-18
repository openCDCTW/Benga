import sqlite3
import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastpCommandline, NcbimakeblastdbCommandline


class SQLiteDatabase:
    def __init__(self, database):
        self._database_path = database

    def connect(self):
        return sqlite3.connect(self._database_path, timeout=30)

    def fetch_scheme(self):
        con = self.connect()
        result = con.execute("select * from scheme").fetchall()
        con.close()
        return result

    def create_table(self):
        con = self.connect()
        con.executescript(
            """
            CREATE TABLE scheme (
                locus_id text PRIMARY KEY,
                dna_seq text
                    );
            CREATE TABLE alleles (
                allele_id text(64) PRIMARY KEY,
                locus_id text not null,
                FOREIGN KEY (locus_id) REFERENCES scheme (locus_id)
                    );
            """
                                       )
        con.close()

    def insert(self, table, values):
        sql = f"insert or replace into {table} values (?, ?)"
        con = self.connect()
        try:
            con.executemany(sql, values)
        except sqlite3.OperationalError:
            pass
        finally:
            con.commit()
            con.close()

    def search(self, query):
        sql = f"select * from alleles where allele_id in ({','.join('?' * len(query))})"
        con = self.connect()
        result = con.execute(sql, query).fetchall()
        con.close()
        return result


def syscall(cmd, stdout=False, stderr=False):
    shell = True if isinstance(cmd, str) else False
    if stdout is False:
        stdout_arg = subprocess.DEVNULL
    elif stdout is None:
        stdout_arg = stderr
    else:
        stdout_arg = subprocess.PIPE
    if stderr is False:
        stderr_arg = subprocess.DEVNULL
    elif stderr is None:
        stderr_arg = stderr
    else:
        stderr_arg = subprocess.PIPE
    child_process = subprocess.run(
        cmd, stdout=stdout_arg, stderr=stderr_arg, check=True, shell=shell, universal_newlines=True
    )
    return child_process


def make_blast_database(input_file, out, dbtype='prot'):
    cline = NcbimakeblastdbCommandline(input_file=input_file, dbtype=dbtype, out=out)
    cline()


def run_blastp(query, db, threads=2):
    cline = NcbiblastpCommandline(
        query=query,
        db=db,
        evalue=1e-6,
        num_threads=threads,
        outfmt='6 qseqid sseqid pident length qlen slen',
    )
    stdout, stderr = cline()
    return stdout, stderr


def generate_record(record_id, sequence, translate=False):
    if isinstance(sequence, Seq) is False:
        sequence = Seq(sequence)
    if translate:
        sequence = sequence.translate(table=11, cds=False)
    return SeqRecord(sequence, id=record_id)
