import psycopg2


with psycopg2.connect(dbname='Salmonella_enterica', host='localhost', user='a504082002', password='i34984047') as conn:
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