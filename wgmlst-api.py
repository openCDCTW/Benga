# This is a RESTful API for profiling

from flask import Flask
from flask_restful import Api
from src.api import profiling


app = Flask(__name__)
api = Api(app)


api.add_resource(profiling.UploadListAPI, '/api/uploads', endpoint='uploads')
api.add_resource(profiling.UploadAPI, '/api/uploads/<string:id>', endpoint='upload')
api.add_resource(profiling.ProfilingAPI, '/api/profiling/<string:id>', endpoint='profiling')
api.add_resource(profiling.ProfileListAPI, '/api/profiles', endpoint='profiles')
api.add_resource(profiling.ProfileAPI, '/api/profiles/<string:id>', endpoint='profile')
api.add_resource(profiling.DendrogramListAPI, '/api/dendrograms', endpoint='dendrograms')
api.add_resource(profiling.DendrogramAPI, '/api/dendrograms/<string:id>', endpoint='dendrogram')

if __name__ == '__main__':
    app.run(debug=True)
