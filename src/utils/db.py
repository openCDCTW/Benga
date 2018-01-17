import json
import subprocess
import psycopg2
import pandas as pd
from sqlalchemy import create_engine, MetaData, Table, Column, ForeignKey
from sqlalchemy.engine.url import URL
from sqlalchemy.dialects import postgresql

DBCONFIG = {}


def load_database_config(filename="database.config", logger=None):
    global DBCONFIG
    with open(filename, "r") as file:
        DBCONFIG = json.loads(file.read())
    DBCONFIG["drivername"] = "postgresql+psycopg2"
    if "host" not in DBCONFIG.keys():
        DBCONFIG["host"] = "localhost"
    if "port" not in DBCONFIG.keys():
        DBCONFIG["port"] = 5432
    logger.info("Read database configuration successfully from {}!".format(filename))
    logger.info("Database: {}:{}".format(DBCONFIG["host"], DBCONFIG["port"]))
    if "password" in DBCONFIG.keys():
        logger.info("Login database as USER {} with PASSWORD ******".format(DBCONFIG["username"]))
    else:
        raise Exception("Passwrod for {} is not available in {}!".format(DBCONFIG["username"], filename))


def from_sql(query, database=None):
    global DBCONFIG
    if database:
        DBCONFIG["database"] = database
    engine = create_engine(URL(**DBCONFIG))
    with engine.connect() as conn:
        t = pd.read_sql_query(query, con=conn)
    engine.dispose()
    return t


def query_filepath(query, database=None):
    global DBCONFIG
    if database:
        DBCONFIG["database"] = database
    with psycopg2.connect(dbname=DBCONFIG["database"], user=DBCONFIG["username"],
                          password=DBCONFIG["password"], host=DBCONFIG["host"],
                          port=DBCONFIG["port"]) as conn:
        with conn.cursor() as cur:
            cur.execute(query)
            filepath = cur.fetchone()[0].tobytes().decode("utf-8")
    return filepath


def to_sql(sql, args={}, database=None):
    global DBCONFIG
    if database:
        DBCONFIG["database"] = database
    engine = create_engine(URL(**DBCONFIG))
    with engine.connect() as conn:
        conn.execute(sql, **args)
    engine.dispose()


def file_to_sql(sql, args, database=None):
    global DBCONFIG
    if database:
        DBCONFIG["database"] = database
    with psycopg2.connect(dbname=DBCONFIG["database"], user=DBCONFIG["username"],
                          password=DBCONFIG["password"], host=DBCONFIG["host"],
                          port=DBCONFIG["port"]) as conn:
        with conn.cursor() as cur:
            cur.execute(sql, args)
            conn.commit()


def append_to_sql(table, df, database=None):
    global DBCONFIG
    if database:
        DBCONFIG["database"] = database
    engine = create_engine(URL(**DBCONFIG))
    with engine.connect() as conn:
        df.to_sql(table, conn, index=False, chunksize=3000, if_exists="append")
    engine.dispose()


def createdb(dbname):
    subprocess.run(["createdb", dbname])


def create_pgadb_relations(dbname, user=None, passwd=None):
    global DBCONFIG
    DBCONFIG["database"] = dbname
    if user:
        DBCONFIG["username"] = user
    if passwd:
        DBCONFIG["password"] = passwd
    engine = create_engine(URL(**DBCONFIG))
    metadata = MetaData()
    loci = Table("loci", metadata,
                 Column("locus_id", None, ForeignKey("locus_meta.locus_id", ondelete="CASCADE"),
                        primary_key=True),
                 Column("ref_allele", None, ForeignKey("alleles.allele_id", ondelete="CASCADE"),
                        nullable=False),
                 Column("occurrence", postgresql.REAL))
    pairs = Table("pairs", metadata,
                    Column("allele_id", None, ForeignKey("alleles.allele_id", ondelete="CASCADE"),
                           primary_key=True, nullable=False),
                    Column("locus_id", None, ForeignKey("locus_meta.locus_id", ondelete="CASCADE"),
                           primary_key=True))
    alleles = Table("alleles", metadata,
                    Column("allele_id", postgresql.CHAR(64), primary_key=True, nullable=False, unique=True),
                    Column("dna_seq", postgresql.TEXT, nullable=False),
                    Column("peptide_seq", postgresql.TEXT, nullable=False),
                    Column("count", postgresql.SMALLINT, nullable=False))
    locus_meta = Table("locus_meta", metadata,
                       Column("locus_id", postgresql.VARCHAR(50), primary_key=True, nullable=False),
                       Column("num_isolates", postgresql.SMALLINT),
                       Column("num_sequences", postgresql.SMALLINT),
                       Column("description", postgresql.TEXT),
                       Column("is_paralog", postgresql.BOOLEAN))
    metadata.create_all(engine)
    engine.dispose()


def create_profiling_relations(user=None, passwd=None):
    global DBCONFIG
    DBCONFIG["database"] = "profiling"
    if user:
        DBCONFIG["username"] = user
    if passwd:
        DBCONFIG["password"] = passwd
    engine = create_engine(URL(**DBCONFIG))
    metadata = MetaData()
    upload = Table("upload", metadata,
                   Column("seq_id", postgresql.CHAR(64), primary_key=True, nullable=False),
                   Column("batch_id", postgresql.CHAR(32), nullable=False),
                   Column("created", postgresql.TIMESTAMP(timezone=True), nullable=False),
                   Column("filename", postgresql.TEXT, nullable=False),
                   Column("file", postgresql.BYTEA, nullable=False))
    profile = Table("profile", metadata,
                    Column("id", postgresql.CHAR(32), primary_key=True, nullable=False),
                    Column("created", postgresql.TIMESTAMP(timezone=True), nullable=False),
                    Column("file", postgresql.TEXT, nullable=False),
                    Column("occurrence", postgresql.SMALLINT, nullable=False),
                    Column("database", postgresql.TEXT, nullable=False))
    dendrogram = Table("dendrogram", metadata,
                       Column("id", postgresql.CHAR(32), primary_key=True, nullable=False),
                       Column("created", postgresql.TIMESTAMP(timezone=True), nullable=False),
                       Column("png_file", postgresql.TEXT, nullable=False),
                       Column("pdf_file", postgresql.TEXT, nullable=False),
                       Column("svg_file", postgresql.TEXT, nullable=False),
                       Column("newick_file", postgresql.TEXT, nullable=False))
    metadata.create_all(engine)
    engine.dispose()
