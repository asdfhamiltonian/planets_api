import os
from flask import Flask, request
from flask_restful import Resource, Api
from collections import OrderedDict

from planets import *
from luna import *

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

class Moon(Resource):
    def get(self, latitude, longitude):
        #flask can't handle negative number input so had to convert from strings
        latitude = float(latitude)
        longitude = float(longitude)
        location = Luna(latitude, longitude)
        location.now()
        moonList = location.list()
        moonDict = OrderedDict()
        moonDict["day"] = moonList[0]
        moonDict["right ascension"] = moonList[1]
        moonDict["declination"] = moonList[2]
        moonDict["latitude"] = moonList[3]
        moonDict["longitude"] = moonList[4]
        moonDict["azimuth"] = moonList[5]
        moonDict["altitude"] = moonList[6]
        moonDict["phase"] = moonList[7]
        moonDict["mean right ascension"] = moonList[8]
        return moonDict

api.add_resource(RiseSet, '/RiseSet/<string:latitude>/<string:longitude>/<string:timezone>')
api.add_resource(PositionNow, '/PositionNow/<string:latitude>/<string:longitude>')
api.add_resource(Moon, '/Moon/<string:latitude>/<string:longitude>')

if __name__ == '__main__':
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port)
