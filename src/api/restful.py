from flask import Flask, abort
from flask_restful import Api, Resource, reqparse

app = Flask(__name__)
api = Api(app)

TASKS = [
    {
        'id': 1,
        'title': u'Buy groceries',
        'description': u'Milk, Cheese, Pizza, Fruit, Tylenol',
        'done': False
    },
    {
        'id': 2,
        'title': u'Learn Python',
        'description': u'Need to find a good Python tutorial on the web',
        'done': False
    }
]


class TaskListAPI(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('title', type=str, required=True,
                                   help='No task title provided', location='json')
        self.reqparse.add_argument('description', type=str, default="", location='json')
        super(TaskListAPI, self).__init__()

    def get(self):
        return {'task': TASKS}

    def post(self):
        task = self.reqparse.parse_args()
        task['id'] = TASKS[-1]['id'] + 1
        TASKS.append(task)
        return TASKS[task['id']-1], 201


class TaskAPI(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('title', type=str, location='json')
        self.reqparse.add_argument('description', type=str, location='json')
        self.reqparse.add_argument('done', type=bool, location='json')
        super(TaskAPI, self).__init__()

    def get(self, id):
        task = [task for task in TASKS if task['id'] == id]
        if len(task) == 0:
            abort(404)
        return {'task': task[0]}


api.add_resource(TaskListAPI, '/api/tasks', endpoint='tasks')
api.add_resource(TaskAPI, '/api/tasks/<int:id>', endpoint='task')

if __name__ == '__main__':
    app.run(debug=True)
