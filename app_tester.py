from flask import Flask
from flask_restful import Resource, Api, reqparse

from gennifer_api_template import generateInputs, run, parseOutput


app = Flask(__name__)
api = Api(app)

parser = reqparse.RequestParser()
parser.add_argument('test')


class RunAlgo(Resource):
    def post(self):
        args = parser.parse_args()
        print(args)
        inputs = generateInputs()
        res = run()
        output = parseOutput()
        return output, 201

api.add_resource(RunAlgo, '/run')

if __name__ == '__main__':
    app.run(debug=True, port=3000)
