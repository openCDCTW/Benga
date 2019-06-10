import os
from django.conf import settings
from pymongo import MongoClient

NOSQLCONFIG = {}


def load_database_config(logger=None):
    global NOSQLCONFIG
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "benga.settings")
    NOSQLCONFIG["host"] = settings.NOSQLS['mongodb']['HOST']
    NOSQLCONFIG["port"] = settings.NOSQLS['mongodb']['PORT']
    logger.info("Database: {}:{}".format(NOSQLCONFIG["host"], NOSQLCONFIG["port"]))


def get_dbtrack(database):
    global NOSQLCONFIG
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "benga.settings")
    NOSQLCONFIG["host"] = settings.NOSQLS['mongodb']['HOST']
    NOSQLCONFIG["port"] = settings.NOSQLS['mongodb']['PORT']
    client = MongoClient(NOSQLCONFIG["host"], NOSQLCONFIG["port"])
    db = client[database]
    return db.track
