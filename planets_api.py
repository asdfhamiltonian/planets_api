import os
from flask import Flask, request
from flask_restful import Resource, Api
from planets import Planets

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
        now = {}
        now["Venus"] = location.calcVenusNow()
        now["Mercury"] = location.calcMercuryNow()
        now["Mars"] = location.calcMarsNow()
        now["Jupiter"] = location.calcJupiterNow()
        now["Saturn"] = location.calcSaturnNow()
        now["Uranus"] = location.calcUranusNow()
        now["Neptune"] = location.calcNeptuneNow()
        now["Pluto"] = location.calcPlutoThatIsNotAPlanetNow()
        return now

api.add_resource(RiseSet, '/RiseSet/<string:latitude>/<string:longitude>/<string:timezone>')
api.add_resource(PositionNow, '/PositionNow/<string:latitude>/<string:longitude>')

if __name__ == '__main__':
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port)
