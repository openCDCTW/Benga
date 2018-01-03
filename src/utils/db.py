import json
import psycopg2
import pandas as pd
import sqlalchemy as sa
from sqlalchemy import create_engine, MetaData, Table, Column, ForeignKey
from sqlalchemy.dialects import postgresql

HOST = "localhost"
PORT = 5432
DATABASE = ""
USER = ""
PASSWORD = ""


def load_database_config(filename="database.config"):
    global HOST
    global PORT
    global DATABASE
    global USER
    global PASSWORD
    with open(filename, "r") as file:
        config = json.loads(file.read())
    DATABASE = config["database"]
    USER = config["user"]
    PASSWORD = config["password"]
    HOST = config["host"] if "host" in config.keys() else HOST
    PORT = config["port"] if "port" in config.keys() else PORT


def from_sql(query, database=None):
    global HOST
    global PORT
    global DATABASE
    global USER
    global PASSWORD
    if not database:
        database = DATABASE
    prot = "postgresql+psycopg2://{}:{}@{}:{}/{}".format(USER, PASSWORD, HOST, PORT, database)
    engine = sa.create_engine(prot)
    with engine.connect() as conn:
        t = pd.read_sql_query(query, con=conn)
    engine.dispose()
    return t


def sql_query_filepath(query, database=None):
    global HOST
    global PORT
    global DATABASE
    global USER
    global PASSWORD
    if not database:
        database = DATABASE
    try:
        with psycopg2.connect(dbname=database, user=USER, password=PASSWORD,
                              host=HOST, port=PORT) as conn:
            with conn.cursor() as cur:
                cur.execute(query)
                return cur.fetchone()[0].tobytes().decode("utf-8")
    except (Exception, psycopg2.DatabaseError) as error:
        print(error)


def to_sql(sql, args, database=None):
    global HOST
    global PORT
    global DATABASE
    global USER
    global PASSWORD
    if not database:
        database = DATABASE
    try:
        with psycopg2.connect(dbname=database, user=USER, password=PASSWORD,
                              host=HOST, port=PORT) as conn:
            with conn.cursor() as cur:
                cur.execute(sql, args)
                conn.commit()
    except (Exception, psycopg2.DatabaseError) as error:
        print(error)


def append_to_sql(df, database=None):
    global HOST
    global PORT
    global DATABASE
    global USER
    global PASSWORD
    if not database:
        database = DATABASE
    try:
        with psycopg2.connect(dbname=database, user=USER, password=PASSWORD,
                              host=HOST, port=PORT) as conn:
            df.to_sql(conn, if_exists="append")
    except (Exception, psycopg2.DatabaseError) as error:
        print(error)


def create_pgadb(dbname, user=None, passwd=None):
    global HOST
    global PORT
    global USER
    global PASSWORD
    if not user:
        user = USER
    if not passwd:
        passwd = PASSWORD
    prot = "postgresql+psycopg2://{}:{}@{}:{}/{}".format(user, passwd, HOST, PORT, dbname)
    engine = create_engine(prot)
    metadata = MetaData()
    loci = Table("loci", metadata,
                 Column("locus_id", postgresql.VARCHAR(50), primary_key=True, nullable=False),
                 Column("ref_allele", None, ForeignKey("alleles.allele_id", ondelete="CASCADE"),
                        nullable=False),
                 Column("occurrence", postgresql.REAL))
    alleles = Table("alleles", metadata,
                    Column("allele_id", postgresql.VARCHAR(16), primary_key=True, nullable=False, unique=True),
                    Column("locus_id", None, ForeignKey("loci.locus_id", ondelete="CASCADE"),
                           primary_key=True),
                    Column("dna_seq", postgresql.TEXT, nullable=False),
                    Column("peptide_seq", postgresql.TEXT, nullable=False),
                    Column("count", postgresql.SMALLINT, nullable=False))
    locus_meta = Table("locus_meta", metadata,
                       Column("locus_id", None, ForeignKey("loci.locus_id", ondelete="CASCADE")),
                       Column("num_isolates", postgresql.SMALLINT),
                       Column("num_sequences", postgresql.SMALLINT),
                       Column("description", postgresql.TEXT),
                       Column("is_paralog", postgresql.BOOLEAN))
    metadata.create_all(engine)