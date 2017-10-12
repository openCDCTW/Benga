import json

import pandas as pd
import sqlalchemy as sa

DATABASE = ""
USER = ""
PASSWORD = ""


def load_database_config(filename="database.config"):
    global DATABASE
    global USER
    global PASSWORD
    with open(filename, "r") as file:
        config = json.loads(file.read())
    DATABASE = config["database"]
    USER = config["user"]
    PASSWORD = config["password"]


def sql_query(query):
    global DATABASE
    global USER
    global PASSWORD
    prot = "postgresql+psycopg2://{}:{}@localhost:5432/{}".format(USER, PASSWORD, DATABASE)
    engine = sa.create_engine(prot)
    conn = engine.connect()
    t = pd.read_sql_query(query, con=conn)
    conn.close()
    engine.dispose()
    return t