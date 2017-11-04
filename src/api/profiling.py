from flask import abort
from flask_restful import Resource, reqparse
from werkzeug.datastructures import FileStorage
import psycopg2
import datetime
import hashlib
import os
from src.utils import db

INDIR = "input"
OUTDIR = "output"
DB = "profiling"


def create_if_not_exist(path):
    if not os.path.exists(path):
        os.makedirs(path)


def get_batch_id(t):
    m = hashlib.md5()
    m.update(t.encode("ascii"))
    return m.hexdigest()


def get_seq_id(file):
    file.seek(0)
    m = hashlib.sha256()
    m.update(file.read())
    return m.hexdigest()


class UploadListAPI(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('filename', type=str, required=True, location='form')
        self.reqparse.add_argument('file', type=FileStorage, required=True, location='files')
        super(UploadListAPI, self).__init__()

    def get(self):
        sql = "select seq_id, filename from upload;"
        results = db.sql_query(sql, database=DB).to_dict(orient="records")
        return results

    def post(self):
        data = self.reqparse.parse_args()
        data["created"] = str(datetime.datetime.now())
        data['batch_id'] = get_batch_id(data["created"])
        file = data.pop('file', None)
        data['seq_id'] = get_seq_id(file)

        input_dir = os.path.join(INDIR, data['batch_id'])
        create_if_not_exist(input_dir)
        file.save(os.path.join(input_dir, data['filename'] + ".fa"))

        sql = "INSERT INTO upload (seq_id,batch_id,created,filename,file) VALUES(%s,%s,%s,%s,%s);"
        args = (data['seq_id'], data['batch_id'], data["created"],
                data['filename'], psycopg2.Binary(file.read()))
        db.to_sql(sql, args, database="profiling")
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


class ProfileListAPI(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('id', type=str, required=True, location='json')
        self.reqparse.add_argument('database', type=str, required=True, location='json')
        self.reqparse.add_argument('occurrence', type=int, required=True, location='json')
        self.reqparse.add_argument('file', type=FileStorage, required=True, location='files')
        super(ProfileListAPI, self).__init__()

    def get(self):
        sql = "select id, occurrence, database from profile;"
        results = db.sql_query(sql, database=DB).to_dict(orient="records")
        return results

    def post(self):
        data = self.reqparse.parse_args()
        data["created"] = str(datetime.datetime.now())
        file = data.pop('file', None)

        sql = "INSERT INTO profile (id,created,file,occurrence,database) VALUES(%s,%s,%s,%s,%s);"
        data = (data['id'], data["created"], psycopg2.Binary(file.read()),
                data["occurrence"], data["database"])
        db.to_sql(sql, data, database="profiling")
        return data, 201


class ProfileAPI(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('id', type=str, required=True, location='json')
        self.reqparse.add_argument('database', type=str, required=True, location='json')
        self.reqparse.add_argument('occurrence', type=int, required=True, location='json')
        self.reqparse.add_argument('file', type=FileStorage, required=True, location='files')
        super(ProfileAPI, self).__init__()

    def get(self, id):
        sql = None
        if len(id) == 32:
            sql = "select id, created, occurrence, database from profile where id='{}';".format(id)
        else:
            abort(404)

        results = db.sql_query(sql, database=DB).to_dict(orient="records")
        if len(results) != 0:
            return results
        else:
            abort(404)


class DendrogramListAPI(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('id', type=str, required=True, location='json')
        self.reqparse.add_argument('png_file', type=FileStorage, required=True, location='files')
        self.reqparse.add_argument('pdf_file', type=FileStorage, required=True, location='files')
        self.reqparse.add_argument('svg_file', type=FileStorage, required=True, location='files')
        self.reqparse.add_argument('newick_file', type=FileStorage, required=True, location='files')
        super(DendrogramListAPI, self).__init__()

    def get(self):
        sql = "select id from dendrogram;"
        results = db.sql_query(sql, database=DB).to_dict(orient="records")
        return results

    def post(self):
        data = self.reqparse.parse_args()
        data["created"] = str(datetime.datetime.now())
        png_file = data.pop('png_file', None)
        pdf_file = data.pop('pdf_file', None)
        svg_file = data.pop('svg_file', None)
        newick_file = data.pop('newick_file', None)

        sql = "INSERT INTO dendrogram (id,created,png_file,pdf_file,svg_file,newick_file) VALUES(%s,%s,%s,%s,%s);"
        data = (data['id'], data["created"], psycopg2.Binary(png_file.read()),
                psycopg2.Binary(pdf_file.read()), psycopg2.Binary(svg_file.read()),
                psycopg2.Binary(newick_file.read()))
        db.to_sql(sql, data, database="profiling")
        return data, 201


class DendrogramAPI(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('id', type=str, required=True, location='json')
        super(DendrogramAPI, self).__init__()

    def get(self, id):
        sql = None
        if len(id) == 32:
            sql = "select id, created from dendrogram where id='{}';".format(id)
        else:
            abort(404)

        results = db.sql_query(sql, database=DB).to_dict(orient="records")
        if len(results) != 0:
            return results
        else:
            abort(404)
