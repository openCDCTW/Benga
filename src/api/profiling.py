from flask import abort, send_file
from flask_restful import Resource, reqparse
from werkzeug.datastructures import FileStorage
import psycopg2
from datetime import datetime
import hashlib
import os
from multiprocessing import Pool
from src.api import internals
from src.utils import db

INDIR = "input"
OUTDIR = "output"
DB = "profiling"
db.load_database_config()


def create_if_not_exist(path):
    if not os.path.exists(path):
        os.makedirs(path)


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
        results = db.sql_query(sql, database=DB).to_dict(orient="records")
        return results

    def post(self):
        data = self.reqparse.parse_args()
        data["created"] = data["created"].strftime('%Y-%m-%d %H:%M:%S')
        file = data.pop('file', None)
        data['seq_id'] = get_seq_id(file)

        input_dir = os.path.join(INDIR, data['batch_id'])
        create_if_not_exist(input_dir)
        file.save(os.path.join(input_dir, data['filename'] + ".fa"))

        sql = "INSERT INTO upload (seq_id,batch_id,created,filename,file) VALUES(%s,%s,%s,%s,%s);"
        args = (data['seq_id'], data['batch_id'], data["created"], data['filename'],
                psycopg2.Binary(file.read()))
        db.to_sql(sql, args, database=DB)
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

        results = db.sql_query(sql, database=DB).to_dict(orient="records")
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
        sql = None
        if len(id) == 32:
            sql = "select * from upload where batch_id='{}';".format(id)
        else:
            abort(404)

        results = db.sql_query(sql, database=DB).to_dict(orient="records")
        if len(results) != 0:
            pool = Pool(processes=1)
            args = [(id, "Salmonella_5k", 95)]
            result = pool.apply_async(internals.profiling_api, args)
            return {"message": "Profiling dataset {}".format(id)}, 200
        else:
            abort(404)


class ProfileListAPI(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('id', type=str, required=True, location='form')
        self.reqparse.add_argument('database', type=str, required=True, location='form')
        self.reqparse.add_argument('occurrence', type=int, required=True, location='form')
        self.reqparse.add_argument('file', type=FileStorage, required=True, location='files')
        super(ProfileListAPI, self).__init__()

    def get(self):
        sql = "select id, occurrence, database from profile;"
        results = db.sql_query(sql, database=DB).to_dict(orient="records")
        return results

    # def post(self):
    #     data = self.reqparse.parse_args()
    #     data["created"] = str(datetime.datetime.now())
    #     file = data.pop('file', None)
    #
    #     sql = "INSERT INTO profile (id,created,file,occurrence,database) VALUES(%s,%s,%s,%s,%s);"
    #     data = (data['id'], data["created"], psycopg2.Binary(file.read()),
    #             data["occurrence"], data["database"])
    #     db.to_sql(sql, data, database=DB)
    #     return data, 201


class ProfileAPI(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('id', type=str, required=True, location='json')
        self.reqparse.add_argument('database', type=str, required=True, location='json')
        self.reqparse.add_argument('occurrence', type=int, required=True, location='json')
        self.reqparse.add_argument('file', type=FileStorage, required=True, location='files')
        super(ProfileAPI, self).__init__()

    def get(self, id):
        if len(id) != 32:
            abort(404)
        sql = "select file from profile where id='{}';".format(id)
        file = db.sql_query_file(sql)
        if not file:
            abort(404)
        return send_file(file, mimetype="text/tab-separated-values")


class DendrogramListAPI(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('id', type=str, required=True, location='form')
        self.reqparse.add_argument('png_file', type=FileStorage, required=True, location='files')
        self.reqparse.add_argument('pdf_file', type=FileStorage, required=True, location='files')
        self.reqparse.add_argument('svg_file', type=FileStorage, required=True, location='files')
        self.reqparse.add_argument('newick_file', type=FileStorage, required=True, location='files')
        super(DendrogramListAPI, self).__init__()

    def get(self):
        sql = "select id from dendrogram;"
        results = db.sql_query(sql, database=DB).to_dict(orient="records")
        return results

    # def post(self):
    #     data = self.reqparse.parse_args()
    #     data["created"] = str(datetime.datetime.now())
    #     png_file = data.pop('png_file', None)
    #     pdf_file = data.pop('pdf_file', None)
    #     svg_file = data.pop('svg_file', None)
    #     newick_file = data.pop('newick_file', None)
    #
    #     sql = "INSERT INTO dendrogram (id,created,png_file,pdf_file,svg_file,newick_file) VALUES(%s,%s,%s,%s,%s);"
    #     data = (data['id'], data["created"], psycopg2.Binary(png_file.read()),
    #             psycopg2.Binary(pdf_file.read()), psycopg2.Binary(svg_file.read()),
    #             psycopg2.Binary(newick_file.read()))
    #     db.to_sql(sql, data, database=DB)
    #     return data, 201


class DendrogramAPI(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('id', type=str, required=True, location='json')
        super(DendrogramAPI, self).__init__()

    def get(self, filetype, id):
        if len(id) != 32:
            abort(404)
        sql = "select {} from dendrogram where id='{}';".format(filetype + "_file", id)
        file = db.sql_query_file(sql)
        if not file:
            abort(404)
        return send_file(file)
