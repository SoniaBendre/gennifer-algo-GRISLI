from flask import Flask
from flask_restful import Resource, Api, reqparse

from gennifer_api import generateInputs, run, parseOutput


app = Flask(__name__)
api = Api(app)

parser = reqparse.RequestParser()
parser.add_argument('exp_dataset_uri')
parser.add_argument('algo')


class RunAlgo(Resource):
    def post(self):
        args = parser.parse_args()
        inputs = generateInputs(args['exp_dataset_uri'])
        res = run(inputs, args['algo'])
        output = parseOutput(*res)
        return output, 201

api.add_resource(RunAlgo, '/run')

if __name__ == '__main__':
    app.run(debug=True)
