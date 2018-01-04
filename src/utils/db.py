import json
import psycopg2
import pandas as pd
from sqlalchemy import create_engine, MetaData, Table, Column, ForeignKey
from sqlalchemy.engine.url import URL
from sqlalchemy.dialects import postgresql

DBCONFIG = {}


def load_database_config(filename="database.config"):
    global DBCONFIG
    with open(filename, "r") as file:
        DBCONFIG = json.loads(file.read())
    DBCONFIG["drivername"] = "postgresql+psycopg2"
    if "host" not in DBCONFIG.keys():
        DBCONFIG["host"] = "localhost"
    if "port" not in DBCONFIG.keys():
        DBCONFIG["port"] = 5432


def from_sql(query, database=None):
    global DBCONFIG
    if database:
        DBCONFIG["database"] = database
    engine = create_engine(URL(DBCONFIG))
    with engine.connect() as conn:
        t = pd.read_sql_query(query, con=conn)
    engine.dispose()
    return t


def sql_query_filepath(query, database=None):
    global DBCONFIG
    if database:
        DBCONFIG["database"] = database
    try:
        with psycopg2.connect(dbname=DBCONFIG["database"], user=DBCONFIG["username"],
                              password=DBCONFIG["password"], host=DBCONFIG["host"],
                              port=DBCONFIG["port"]) as conn:
            with conn.cursor() as cur:
                cur.execute(query)
                return cur.fetchone()[0].tobytes().decode("utf-8")
    except (Exception, psycopg2.DatabaseError) as error:
        print(error)


def to_sql(sql, args, database=None):
    global DBCONFIG
    if database:
        DBCONFIG["database"] = database
    try:
        engine = create_engine(URL(DBCONFIG))
        with engine.connect() as conn:
            conn.execute(sql, **args)
        engine.dispose()
    except Exception as error:
        print(error)


def file_to_sql(sql, args, database=None):
    global DBCONFIG
    if database:
        DBCONFIG["database"] = database
    try:
        with psycopg2.connect(dbname=DBCONFIG["database"], user=DBCONFIG["username"],
                              password=DBCONFIG["password"], host=DBCONFIG["host"],
                              port=DBCONFIG["port"]) as conn:
            with conn.cursor() as cur:
                cur.execute(sql, args)
                conn.commit()
    except (Exception, psycopg2.DatabaseError) as error:
        print(error)


def append_to_sql(table, df, database=None):
    global DBCONFIG
    if database:
        DBCONFIG["database"] = database
    try:
        engine = create_engine(URL(DBCONFIG))
        with engine.connect() as conn:
            df.to_sql(table, conn, index=False, chunksize=3000, if_exists="append")
    except (Exception, psycopg2.DatabaseError) as error:
        print(error)
    finally:
        engine.dispose()


def create_pgadb(dbname, user=None, passwd=None):
    global DBCONFIG
    DBCONFIG["database"] = dbname
    if user:
        DBCONFIG["username"] = user
    if passwd:
        DBCONFIG["password"] = passwd
    engine = create_engine(URL(DBCONFIG))
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
    engine.dispose()
