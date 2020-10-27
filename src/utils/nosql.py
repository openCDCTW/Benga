import os
from django.conf import settings
from pymongo import MongoClient

NOSQLCONFIG = {}


class NoSQL:
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "benga.settings")

    def __init__(self, database_name):
        global NOSQLCONFIG
        os.environ.setdefault("DJANGO_SETTINGS_MODULE", "benga.settings")
        NOSQLCONFIG["host"] = settings.NOSQLS['mongodb']['HOST']
        NOSQLCONFIG["port"] = int(settings.NOSQLS['mongodb']['PORT'])
        self.client = MongoClient(NOSQLCONFIG["host"], NOSQLCONFIG["port"])
        self.database = self.client[database_name]

    def disconnect(self):
        self.client.close()

    def profilesIterator(self):
        collection = self.database['Sample_attribute']
        return collection.find({}, {'profile': 1, 'NCBIAccession': 1})

    @property
    def fields(self):
        return self.database.Description.find_one({'fields': {'$exists': True}})['fields']

    @property
    def core_genome(self):
        return self.database.Description.find_one({'core_genome': {'$exists': True}})['core_genome']

    def fetch_attrib(self, filter_, skip):
        return self.database.Sample_attribute.find(filter_, skip)

