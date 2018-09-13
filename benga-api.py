from flask import Flask
from flask_restful import Api

from benga.src.api import profiling

app = Flask(__name__)
api = Api(app)


api.add_resource(profiling.ProfilingAPI, '/api/profiling/<string:id>', endpoint='profiling')

if __name__ == '__main__':
    app.run(debug=True)
