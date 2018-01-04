import json
import psycopg2
import pandas as pd
import sqlalchemy as sa

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


def sql_query(query, database=None):
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
