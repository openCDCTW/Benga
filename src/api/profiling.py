from flask import abort, send_file
from flask_restful import Resource, reqparse
from werkzeug.datastructures import FileStorage
import psycopg2
from datetime import datetime
import hashlib
import os
from threading import Thread
from src.api import internals
from src.utils import db
from src.utils.files import create_if_not_exist
from src.models import logs

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


class UploadListAPI(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('batch_id', type=str, required=True, location='form')
        self.reqparse.add_argument('created', type=lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S'),
                                   required=True, location='form')
        self.reqparse.add_argument('filename', type=str, required=True, location='form')
        self.reqparse.add_argument('file', type=FileStorage, required=True, location='files')
        super(UploadListAPI, self).__init__()

    def get(self):
        sql = "select seq_id, filename from upload;"
        results = db.from_sql(sql, database=DB).to_dict(orient="records")
        return results

    def post(self):
        data = self.reqparse.parse_args()
        input_dir = os.path.join(INDIR, data['batch_id'])
        create_if_not_exist(input_dir)
        data['file'].save(os.path.join(input_dir, data['filename'] + ".fa"))
        data["created"] = data["created"].strftime('%Y-%m-%d %H:%M:%S')
        data['seq_id'] = get_seq_id(data['file'])

        sql = "INSERT INTO upload (seq_id,batch_id,created,filename,file) VALUES(%s,%s,%s,%s,%s);"
        args = (data['seq_id'], data['batch_id'], data["created"], data['filename'],
                psycopg2.Binary(data['file'].read()))
        db.file_to_sql(sql, args, database=DB)
        data.pop('file', None)
        return data, 201


class UploadAPI(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('seq_id', type=str, location='json')
        self.reqparse.add_argument('batch_id', type=str, location='json')
        super(UploadAPI, self).__init__()

    def get(self, id):
        sql = None
        if len(id) == 32:
            sql = "select seq_id, batch_id, filename from upload where batch_id='{}';".format(id)
        elif len(id) == 64:
            sql = "select seq_id, batch_id, filename from upload where seq_id='{}';".format(id)
        else:
            abort(404)

        results = db.from_sql(sql, database=DB).to_dict(orient="records")
        if len(results) != 0:
            return results
        else:
            abort(404)


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


class ProfileListAPI(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('id', type=str, required=True, location='form')
        self.reqparse.add_argument('database', type=str, required=True, location='form')
        self.reqparse.add_argument('occurrence', type=int, required=True, location='form')
        self.reqparse.add_argument('file', type=str, required=True, location='files')
        super(ProfileListAPI, self).__init__()

    def get(self):
        sql = "select id, occurrence, database from profile;"
        results = db.from_sql(sql, database=DB).to_dict(orient="records")
        return results


class ProfileAPI(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('id', type=str, required=True, location='json')
        self.reqparse.add_argument('database', type=str, required=True, location='json')
        self.reqparse.add_argument('occurrence', type=int, required=True, location='json')
        self.reqparse.add_argument('file', type=str, required=True, location='files')
        super(ProfileAPI, self).__init__()

    def get(self, id):
        if len(id) != 32:
            abort(404)
        sql = "select file from profile where id='{}';".format(id)
        filepath = db.query_filepath(sql, database="profiling")
        if not filepath:
            abort(404)
        return send_file(filepath, mimetype="text/tab-separated-values")


class DendrogramListAPI(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('id', type=str, required=True, location='form')
        self.reqparse.add_argument('png_file', type=str, required=True, location='files')
        self.reqparse.add_argument('pdf_file', type=str, required=True, location='files')
        self.reqparse.add_argument('svg_file', type=str, required=True, location='files')
        self.reqparse.add_argument('newick_file', type=str, required=True, location='files')
        super(DendrogramListAPI, self).__init__()

    def get(self):
        sql = "select id from dendrogram;"
        results = db.from_sql(sql, database=DB).to_dict(orient="records")
        return results


class DendrogramAPI(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('id', type=str, required=True, location='json')
        super(DendrogramAPI, self).__init__()

    def get(self, filetype, id):
        if len(id) != 32:
            abort(404)
        sql = "select {} from dendrogram where id='{}';".format(filetype + "_file", id)
        filepath = db.query_filepath(sql, database="profiling")
        if not filepath:
            abort(404)
        return send_file(filepath)
