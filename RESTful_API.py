from flask import Flask
from flask_restful import Api
from src.api import profiling


app = Flask(__name__)
api = Api(app)


api.add_resource(profiling.UploadListAPI, '/api/uploads', endpoint='uploads')
api.add_resource(profiling.UploadAPI, '/api/uploads/<string:id>', endpoint='upload')

if __name__ == '__main__':
    app.run(debug=True)
