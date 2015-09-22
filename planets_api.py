import os
from flask import Flask, request
from flask_restful import Resource, Api
from planets import Planets
from collections import OrderedDict

app = Flask(__name__)
api = Api(app)

class RiseSet(Resource):
    def get(self, latitude, longitude, timezone):
        #flask can't handle negative number input so had to convert from strings
        latitude = float(latitude)
        longitude = float(longitude)
        timezone = float(timezone)
        location = Planets(latitude, longitude)
        riseSetDict = location.dictRiseSet(-7)
        return riseSetDict

class PositionNow(Resource):
    def get(self, latitude, longitude):
        #flask can't handle negative number input so had to convert from strings
        latitude = float(latitude)
        longitude = float(longitude)
        location = Planets(latitude, longitude)
        now = OrderedDict()
        now["mercury"] = location.calcMercuryNow()
        now["venus"] = location.calcVenusNow()
        now["mars"] = location.calcMarsNow()
        now["jupiter"] = location.calcJupiterNow()
        now["saturn"] = location.calcSaturnNow()
        now["uranus"] = location.calcUranusNow()
        now["neptune"] = location.calcNeptuneNow()
        now["pluto"] = location.calcPlutoThatIsNotAPlanetNow()
        return now

api.add_resource(RiseSet, '/RiseSet/<string:latitude>/<string:longitude>/<string:timezone>')
api.add_resource(PositionNow, '/PositionNow/<string:latitude>/<string:longitude>')

if __name__ == '__main__':
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port)
