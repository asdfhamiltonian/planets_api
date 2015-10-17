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
        # flask can't handle negative numbers so converted from string
        latitude = float(latitude)
        longitude = float(longitude)
        timezone = float(timezone)
        location = Planets(latitude, longitude)
        riseSetDict = location.dictRiseSet(-7)
        return riseSetDict


class PositionNow(Resource):
    def get(self, latitude, longitude):
        # flask can't handle negative numbers so converted from string
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
    def get(self, latitude, longitude, timezone):
        # flask can't handle negative numbers so converted from string
        latitude = float(latitude)
        longitude = float(longitude)
        timezone = float(timezone)
        location = Luna(latitude, longitude)
        moonList = location.now()
        moonDict = OrderedDict()
        moonDict["azimuth"] = moonList[0]
        moonDict["altitude"] = moonList[1]
        moonDict["right ascension"] = moonList[2]
        moonDict["declination"] = moonList[3]
        moonDict["phase"] = moonList[4]
        moonDict["topocentric right ascension"] = moonList[5]
        moonDict["topocentric declination"] = moonList[6]
        moonDict["mpar"] = moonList[7]
        moonDict["msd"] = moonList[8]
        moonDict["longitude degrees"] = moonList[9]
        moonDict["latitude degrees"] = moonList[10]
        moonDict["gmsto"] = moonList[11]
        moonDict["epoch day"] = moonList[12]
        moonDict["riseset1"] = location.risesetlist(1, timezone)[0:3]
        moonDict["riseset2"] = location.risesetlist(2, timezone)[0:3]
        moonDict["riseset3"] = location.risesetlist(3, timezone)[0:3]
        return moonDict

api.add_resource(RiseSet, '/RiseSet/<string:latitude>/<string:longitude>/<string:timezone>')
api.add_resource(PositionNow, '/PositionNow/<string:latitude>/<string:longitude>')
api.add_resource(Moon, '/Moon/<string:latitude>/<string:longitude>/<string:timezone>')

if __name__ == '__main__':
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port)
