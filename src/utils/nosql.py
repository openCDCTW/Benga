import os
from django.conf import settings
from pymongo import MongoClient

NOSQLCONFIG = {}


class NoSQL:
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "benga.settings")

    def __init__(self, logger=None):
        self.config = dict()
        self.config["host"] = settings.NOSQLS['mongodb']['HOST']
        self.config["port"] = int(settings.NOSQLS['mongodb']['PORT'])
        logger.info("Database: {}:{}".format(self.config["host"], self.config["port"]))

    def connect(self, database):
        client = MongoClient(self.config["host"], self.config["port"])
        db = client[database]
        return db['track_data'], db['metadata_cols']


def load_database_config(logger=None):
    global NOSQLCONFIG
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "benga.settings")
    NOSQLCONFIG["host"] = settings.NOSQLS['mongodb']['HOST']
    NOSQLCONFIG["port"] = int(settings.NOSQLS['mongodb']['PORT'])
    logger.info("Database: {}:{}".format(NOSQLCONFIG["host"], NOSQLCONFIG["port"]))


def get_dbtrack(database):
    global NOSQLCONFIG
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "benga.settings")
    NOSQLCONFIG["host"] = settings.NOSQLS['mongodb']['HOST']
    NOSQLCONFIG["port"] = int(settings.NOSQLS['mongodb']['PORT'])
    client = MongoClient(NOSQLCONFIG["host"], NOSQLCONFIG["port"])
    db = client[database]
    return db['track_data'], db['metadata_cols']
