import argparse
import psycopg2
import os


def main():
    args = parse_args()

    dbname = args.database
    user = args.user
    passwd = args.passwd

    try:
        os.system("createdb {}".format(dbname))
        os.system("createdb profiling")
    except:
        print("Failed creating database: {}".format(dbname))

    try:
        create_relations(dbname, user, passwd)
        create_profiling_relations(user, passwd)
    except:
        print("Failed creating relations")


def parse_args():
    arg_parser = argparse.ArgumentParser()

    arg_parser.add_argument(
        "-u", "--user",
        required=True,
        help="Database user. (necessary)"
    )

    arg_parser.add_argument(
        "-p", "--passwd",
        required=True,
        help="Login password to the user. (necessary)"
    )

    arg_parser.add_argument(
        "-d", "--database",
        required=True,
        help="Database name to create. (necessary)"
    )
    return arg_parser.parse_args()


def create_relations(dbname, user, passwd):
    with psycopg2.connect(dbname=dbname, host='localhost', user=user, password=passwd) as conn:
        with conn.cursor() as cur:
            cur.execute("""CREATE TABLE loci (
                locus_id text NOT NULL PRIMARY KEY,
                ref_allele char(16) NOT NULL,
                occurence real
            );""")

            cur.execute("""CREATE TABLE alleles (
                allele_id char(16) NOT NULL PRIMARY KEY,
                locus_id text references loci(locus_id),
                allele_int integer NOT NULL
            );""")

            cur.execute("ALTER TABLE loci ADD FOREIGN KEY (ref_allele) REFERENCES alleles(allele_id);")

            cur.execute("""CREATE TABLE locus_meta (
                locus_id text NOT NULL references loci(locus_id),
                num_isolates smallserial,
                num_sequences smallserial,
                description text,
                is_paralog boolean
            );""")

            cur.execute("""CREATE TABLE allele_count (
                allele_id char(16) NOT NULL references alleles(allele_id),
                count smallserial
            );""")

            cur.execute("""CREATE TABLE sequences (
                sequence_id char(16) NOT NULL references alleles(allele_id),
                dna_seq text NOT NULL,
                peptide_seq text NOT NULL
            );""")


def create_profiling_relations(user, passwd):
    with psycopg2.connect(dbname="profiling", host='localhost', user=user, password=passwd) as conn:
        with conn.cursor() as cur:
            cur.execute("""CREATE TABLE upload (
                seq_id char(64) NOT NULL PRIMARY KEY,
                batch_id char(32) NOT NULL,
                created timestamp with time zone NOT NULL,
                filename text NOT NULL,
                file bytea NOT NULL
            );""")

            cur.execute("""CREATE TABLE profile (
                id char(32) NOT NULL PRIMARY KEY,
                created timestamp with time zone NOT NULL,
                file text NOT NULL,
                occurrence smallint NOT NULL,
                database text NOT NULL
            );""")

            cur.execute("""CREATE TABLE dendrogram (
                id char(32) NOT NULL PRIMARY KEY,
                created timestamp with time zone NOT NULL,
                png_file text NOT NULL,
                pdf_file text NOT NULL,
                svg_file text NOT NULL,
                newick_file text NOT NULL
            );""")


if __name__ == "__main__":
    main()
