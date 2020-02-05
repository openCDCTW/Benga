import os
import subprocess
import pandas as pd
from sqlalchemy import create_engine, MetaData, Table, Column, ForeignKey
from sqlalchemy.engine.url import URL
from sqlalchemy.dialects import postgresql
from sqlalchemy_utils.functions import drop_database
from django.conf import settings

DBCONFIG = {}
METADATA = MetaData()
LOCI = Table("loci", METADATA,
             Column("locus_id", None, ForeignKey("locus_meta.locus_id", ondelete="CASCADE"),
                    primary_key=True),
             Column("ref_allele", None, ForeignKey("alleles.allele_id", ondelete="CASCADE"),
                    nullable=False),
             Column("occurrence", postgresql.REAL))
PAIRS = Table("pairs", METADATA,
              Column("allele_id", None, ForeignKey("alleles.allele_id", ondelete="CASCADE"),
                     primary_key=True, nullable=False),
              Column("locus_id", None, ForeignKey("locus_meta.locus_id", ondelete="CASCADE"),
                     primary_key=True))
ALLELES = Table("alleles", METADATA,
                Column("allele_id", postgresql.CHAR(64), primary_key=True, nullable=False, unique=True),
                Column("dna_seq", postgresql.TEXT, nullable=False),
                Column("peptide_seq", postgresql.TEXT, nullable=False),
                Column("count", postgresql.SMALLINT, nullable=False))
LOCUS_META = Table("locus_meta", METADATA,
                   Column("locus_id", postgresql.VARCHAR(50), primary_key=True, nullable=False),
                   Column("num_isolates", postgresql.SMALLINT),
                   Column("num_sequences", postgresql.SMALLINT),
                   Column("description", postgresql.TEXT),
                   Column("is_paralog", postgresql.BOOLEAN))


def load_database_config(logger=None):
    global DBCONFIG
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "benga.settings")
    DBCONFIG["drivername"] = "postgresql+psycopg2"
    DBCONFIG["host"] = settings.DATABASES['default']['HOST']
    DBCONFIG["port"] = settings.DATABASES['default']['PORT']
    DBCONFIG["username"] = settings.DATABASES['default']['USER']
    DBCONFIG["password"] = settings.DATABASES['default']['PASSWORD']
    logger.info("Database: {}:{}".format(DBCONFIG["host"], DBCONFIG["port"]))
    logger.info("Login database as USER {} with PASSWORD ******".format(DBCONFIG["username"]))


def from_sql(query, database=None):
    global DBCONFIG
    if database:
        DBCONFIG["database"] = database
    engine = create_engine(URL(**DBCONFIG))
    with engine.connect() as conn:
        t = pd.read_sql_query(query, con=conn)
    engine.dispose()
    return t


def to_sql(sql, args={}, database=None):
    global DBCONFIG
    if database:
        DBCONFIG["database"] = database
    engine = create_engine(URL(**DBCONFIG))
    with engine.connect() as conn:
        conn.execute(sql, **args)
    engine.dispose()


def table_to_sql(table, df, database=None, append=True):
    global DBCONFIG
    if database:
        DBCONFIG["database"] = database
    if_exists = "append" if append else "fail"
    engine = create_engine(URL(**DBCONFIG))
    with engine.connect() as conn:
        df.to_sql(table, conn, index=False, chunksize=3000, if_exists=if_exists)
    engine.dispose()


def createdb(dbname):
    subprocess.run(["createdb", dbname])


def create_pgadb_relations(dbname):
    global DBCONFIG
    DBCONFIG["database"] = dbname
    engine = create_engine(URL(**DBCONFIG))
    METADATA.create_all(engine)
    engine.dispose()


def dropdb(dbname):
    global DBCONFIG
    DBCONFIG["database"] = dbname
    engine = create_engine(URL(**DBCONFIG))
    drop_database(engine.url)
