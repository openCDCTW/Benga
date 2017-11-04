from flask import Flask, abort
from flask_restful import Api, Resource, reqparse
from werkzeug.datastructures import FileStorage
import datetime
import hashlib
import os
from src.utils import db

INDIR = "input"
OUTDIR = "output"
DB = "profiling"
UPLOADS = []

app = Flask(__name__)
api = Api(app)


def create_if_not_exist(path):
    if not os.path.exists(path):
        os.makedirs(path)


def get_batch_id(t):
    m = hashlib.md5()
    m.update(t.encode())
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
        sql = "select filename from upload;"
        results = db.sql_query(sql, database=DB)
        return results
        # return UPLOADS

    def post(self):
        upload = self.reqparse.parse_args()
        upload["created"] = str(datetime.datetime.now())
        upload['batch_id'] = get_batch_id(upload["created"])
        input_dir = os.path.join(INDIR, upload['batch_id'])
        create_if_not_exist(input_dir)
        file = upload.pop('file', None)
        file.save(os.path.join(input_dir, upload['filename'] + ".fa"))
        upload['seq_id'] = get_seq_id(file)
        UPLOADS.append(upload)
        return upload, 201


class UploadAPI(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('seq_id', type=str, location='json')
        self.reqparse.add_argument('batch_id', type=str, location='json')
        super(UploadAPI, self).__init__()

    def get(self, id):
        item = None
        if len(id) == 32:
            item = 'batch_id'
        elif len(id) == 64:
            item = 'seq_id'
        else:
            abort(404)
        for each in UPLOADS:
            if each[item] == id:
                return each
        else:
            abort(404)


api.add_resource(UploadListAPI, '/api/uploads', endpoint='uploads')
api.add_resource(UploadAPI, '/api/uploads/<int:id>', endpoint='upload')

if __name__ == '__main__':
    app.run(debug=True)
