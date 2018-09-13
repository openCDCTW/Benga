import hashlib
from threading import Thread
from flask import abort
from flask_restful import Resource, reqparse
from benga.src.api import internals
from benga.src.utils import db, logs
from benga.src.utils.files import create_if_not_exist

INDIR = "input"
OUTDIR = "output"
DB = "profiling"

create_if_not_exist(OUTDIR)
lf = logs.LoggerFactory()
lf.addConsoleHandler()
logger = lf.create()
db.load_database_config(logger=logger)


def get_seq_id(file):
    file.seek(0)
    m = hashlib.sha256()
    m.update(file.read())
    return m.hexdigest()

class ProfilingAPI(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('id', type=str, required=True, location='json')
        self.reqparse.add_argument('database', type=str, required=True, location='json')
        self.reqparse.add_argument('occurrence', type=int, required=True, location='json')
        super(ProfilingAPI, self).__init__()

    def get(self, id):
        if len(id) != 32:
            abort(404)
        sql = "select * from upload where batch_id='{}';".format(id)
        results = db.from_sql(sql, database=DB).to_dict(orient="records")
        if len(results) == 0:
            abort(404)
        Thread(target=internals.profiling_api, args=(id, "Salmonella_5k", 95), daemon=True).start()
        return {"message": "Profiling dataset {}".format(id)}, 200

