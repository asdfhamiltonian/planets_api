"""
Moon 4 arcminute calculation
Adapted from original code in QBASIC by Keith Burnett
available through this website: http://www.stargazing.net/kepler/moon2.html
Futher adatped with algorithms/equations from page "How to Compute Planetary Positions"
by Paul Schlyter http://www.stjarnhimlen.se/comp/ppcomp.html#20
"""

from math import *
import datetime

def degrees(rads):
        return (rads/pi)*180

def radians(degrees):
    return (degrees/180 * pi)

def fnday(y,m,d,h):
    a = 367 * y - 7 * (y + (m + 9) // 12) // 4 + 275 * m // 9 + d - 730530 + h / 24
    return a

def fnipart(x):
    return copysign(1,x) * floor(abs(x))

def fnrange(angle):
    newangle = angle%(2*pi)
    if newangle < 0:
        newangle = newangle + 2*pi
    return newangle

def fnatan2(y, x):
    a = atan2(y, x)
    if a < 0:
        a = a + 2 * pi
    return a

def timeadjust(x, utc):
    x = x + utc
    while x < 0:
        x = x + 24
    if x > 24:
        x = x - 24
    return x

def decimaltominsecs(x):
    x = x % 1
    mins = floor(x * 60)
    x = (x * 60) % 1
    secs = floor(x * 60)
    return (mins, secs)

def longituderange(phi):
    return (phi + 180)%360 - 180

def latituderange(theta):
    return abs(((theta - 90)%360) - 180) - 90

class Luna(object):
    localLat = 0
    localLong = 0
    d = 0
    ra = 0
    dec = 0
    lstring = ""
    lonstring = ""
    az = 0
    alttopoc = 0
    phasestr = ""
    mra = 0
    mdecl = 0
    mpar = 0
    gmsto = 0
    msd = 0
    
    
    def __init__(self, llat = 45, llon = -122):
        self.localLat = llat
        self.localLong = llon
        
    def calculate(self, y, m, d, h, mins, second):
        """
        y = datetime.datetime.utcnow().year
        m = datetime.datetime.utcnow().month
        d = datetime.datetime.utcnow().day
        h = datetime.datetime.utcnow().hour
        mins = datetime.datetime.utcnow().minute
        second = datetime.datetime.utcnow().second
        """
        
        h = h + mins / 60 + second / 3600
        d = fnday(y,m,d,h)
        
        #These are the moon elements:

        Nm = fnrange(radians(125.1228 - .0529538083 * d))
        im = radians(5.1454)
        wm = fnrange(radians(318.0634 + 0.1643573223 * d))
        am = 60.2666 #mean radius in earth radii?
        ecm = 0.0549
        Mm = fnrange(radians(115.3654 + 13.0649929509 * d))

        #Sun elements are next:
        Ns = 0
        isun = 0
        ws = fnrange(radians(282.9404 + 4.70935E-05 * d))
        asun = 1.0 #Astronomical units
        ecs = 0.016709 - 1.151E-09 * d
        Ms = fnrange(radians(356.047 + 0.9856002585 * d))

        #sun position (for azimuth info)
        Es = Ms + ecs * sin(Ms) * (1.0 + ecs * cos(Ms))
        xv = cos(Es) - ecs
        yv = sqrt(1 - ecs * ecs) * sin(Es)
        vs = atan2(yv, xv)
        rs = sqrt(xv*xv + yv*yv )
        lonsun = ws + vs
        xs = rs * cos(lonsun)
        ys = rs * sin(lonsun)
        #since the SUn always is in the ecliptic plane, zs is of course zero??
        xes = xs
        yes = ys * cos(ecs)
        zes = ys * sin(ecs)
        ra = atan2(yes, xes)
        dec = atan2(zes, sqrt(xes*xes + yes*yes))

        """
        'position of Sun
        'Es = Ms + Es * SIN(Ms) * (1! + ecs * COS(Ms))
        'xv = COS(Es) - ecs
        'yv = SQR(1! - ecs * ecs) * SIN(Es)
        'vs = FNatn2(yv, xv)
        'rs = SQR(xv * xv + yv * yv)
        'lonsun = vs + ws
        'xs = rs * COS(lonsun)
        'ys = rs * SIN(lonsun)
        'xe = xs
        'ye = ys * COS(ecl)
        ''ze = ys * SIN(ecl)
        'ras = FNatn2(ye, xe)
        'decs = FNatn2(ze, sqr(xe * xe + ye * ye))
        """

        #Position of Moon
        Em = Mm + ecm * sin(Mm) * (1 + ecm * cos(Mm))
        xv = am * (cos(Em) - ecm)
        yv = am * (sqrt(1 - ecm * ecm) * sin(Em))
        vm = fnatan2(yv, xv)
        rm = sqrt(xv * xv + yv * yv)
        xh = rm * (cos(Nm) * cos(vm + wm) - sin(Nm) * sin(vm + wm) * cos(im))
        yh = rm * (sin(Nm) * cos(vm + wm) + cos(Nm) * sin(vm + wm) * cos(im))
        zh = rm * (sin(vm + wm) * sin(im))

        #moon's geometric longitude and latitude:
        lon = fnatan2(yh, xh)
        lat = fnatan2(zh, sqrt(xh * xh + yh * yh))

        """
        perturbations
        first needed to calc the args below, which are in radians
        MS, Mm - Mean anomaly of the Sun and Moon
        Nm - Longitude of the Moon's node
        ws, wm - Argument of the perihelion for the Sun and the Moon
        """

        Ls = Ms + ws # Mean longitude of the Sun (Ns = 0)
        Lm = Mm + wm + Nm #Mean longitude of the Moon
        dm = Lm - Ls # Mean elongatio22n of the Moon
        F = Lm - Nm #Argument of latitude for the Moon

        #Then add the following terms to the longitude
        # Note amplitudes are in degrees, convert at end
        dlon = -1.274 * sin(Mm - 2 * dm) #the Evection?
        dlon = dlon + 0.658 * sin(2 * dm) # the Variation?
        dlon = dlon - 0.186 * sin(Ms) # the Yearly Equation?
        dlon = dlon - 0.059 * sin(2 * Mm - 2 * dm)
        dlon = dlon - 0.057 * sin(Mm - 2*dm + Ms)
        dlon = dlon + 0.053 * sin(Mm + 2 * dm)
        dlon = dlon + 0.046 * sin(2 * dm - Ms)
        dlon = dlon + 0.041 * sin(Mm - Ms)
        dlon = dlon - 0.035 * sin(dm)   #The Parallactic Equation??
        dlon = dlon - 0.031 * sin(Mm + Ms)
        dlon = dlon - 0.015 * sin(2 * F - 2 * dm)
        dlon = dlon + 0.011 * sin(Mm - 4 * dm)
        lon = radians(dlon) + lon
        
        
        #latitude terms
        dlat = -0.173 * sin(F - 2 * dm)
        dlat = dlat - 0.055 * sin(Mm - F - 2 * dm)
        dlat = dlat - 0.046 * sin(Mm + F - 2 * dm)
        dlat = dlat + 0.033 * sin(F + 2 * dm)
        dlat = dlat + 0.017 * sin(2 * Mm + F)
        lat = radians(dlat) + lat

        """
        Need to correct so latitude is between pi/2 and -pi/2 (+90 or -90 degrees)
        The above formula had been giving latitudes like "355.679..."
        There's some math involved that makes sense if you plot desired output vs. input
        """
        
        while lat < 0:
                lat = lat + 2*pi
        
        if lat > pi/2 and lat <= 3*pi/2:
                lat = pi - lat
        elif lat > 3*pi/2 and lat <= 5*pi/2:
                lat = lat - 2*pi
        
        #distance terms earth radii
        rm = rm - 0.58 * cos(Mm - 2 * dm)
        rm = rm - 0.46 * cos(2 * dm)

        #Next find the cartesian coordinates of the geocentric lunar position

        xg = rm * cos(lon) * cos(lat)
        yg = rm * sin(lon) * cos(lat)
        zg = rm * sin(lat)
        #rotate the equatorial coords
        #obliquity of ecliptic of date
        ecl = radians(23.4393 - 3.563E-07 * d)
        xe = xg
        ye = yg * cos(ecl) - zg * sin(ecl)
        ze = yg * sin(ecl) + zg * cos(ecl)

        #geocentric RA and Dec
        ra = fnatan2(ye, xe)
        dec = fnatan2(ze, sqrt(xe*xe + ye*ye))
        lstring = "Lat: " + str(latituderange(degrees(lat)))
        lonstring = "Lon: " + str(longituderange(degrees(lon)))

        """
        next sidereal time:
        print(Ls)
        print(vs+ws)
        localLong = float(input("longitude? "))
        localLat = float(input("latitude? "))
        """
        
        lat1 = radians(self.localLat)
        gmsto = 12 * (Ls + pi)/(pi)
        gmst = gmsto + h
        lst = gmst + self.localLong/15
        lst = lst%24
        ha = lst - degrees(ra)/15
        while (ha > 24) or (ha < -24):
            if ha > 12:
                ha = ha - 24
            elif ha < -12:
                ha = ha + 24
        ha = radians(ha*15)
        x = cos(ha) * cos(dec)
        y = sin(ha) * cos(dec)
        z = sin(dec)

        xhor = x * sin(lat1) - z * cos(lat1)
        yhor = y
        zhor = x * cos(lat1) + z * sin(lat1)

        az = atan2( yhor, xhor) + pi
        alt = asin(zhor)

        #print("Azimuth is {} degrees".format(str(degrees(az))))
        #print("Alt is {} degrees".format(str(degrees(alt))))

        #now to correct for topocentricity?
        #see http://www.stjarnhimlen.se/comp/ppcomp.html#12b for more info about this
        mpar = asin(1 / rm)
        msd = asin(0.2725 / rm)
        alttopoc = alt - mpar * cos(alt) 

        #calculating "geocentric lattitude" accounting for flattening of earth:
        gclat = lat1 - radians(0.1924) * sin(2 * lat1)
        rho = 0.99833 + 0.00167 * cos(2 * lat1)
        g = atan( tan(gclat) / cos(ha)) #g is the auxiliary angle?
        topRA = ra - mpar * rho * cos(gclat) * sin(ha) / cos(dec)
        if abs(abs(dec) - pi/2) < 0.001 or abs(gclat) < 0.001:
            topDecl = dec - mpar * rho * sin(-dec) * cos(ha)
        else:
            topDecl = dec - mpar * rho * sin(gclat) * sin(g - dec) / sin(g)

        #moon phase calcs
        slon = lonsun
        elong = acos( cos(slon - lon) * cos(lat) )
        FV = pi - elong
        phase = (1 + cos(FV))/2
        phasestr = "The current phase of the moon is {}.".format(str(phase))

        self.d = d
        self.ra = ra
        self.dec = dec 
        self.lstring = lstring
        self.lonstring = lonstring 
        self.az = az
        self.alttopoc = alttopoc 
        self.phasestr = phasestr
        self.mra = topRA
        self.mdecl = topDecl
        self.mpar = mpar
        self.msd = msd
        self.gmsto = gmsto

    def list(self):
        """
        Returns current calculated values for the moon position
        as an array in the following format:
        [day, right_ascension, declination, latitude_info, longitude_info, azimuth,
         altitude, phase_info, mean_right_ascension, topocentric_declination, mpar?,
         msd, gmsto]
        """
        return [self.d, self.ra, self.dec, self.lstring, self.lonstring, self.az,
                self.alttopoc, self.phasestr, self.mra, self.mdecl, self.mpar,
                self.msd, self.gmsto]

    def calclist(self, y, m, d, h, mins, second):
        self.calculate(y, m, d, h, mins, second)
        return self.list()
        
    def now(self):
        y = datetime.datetime.utcnow().year
        m = datetime.datetime.utcnow().month
        d = datetime.datetime.utcnow().day
        h = datetime.datetime.utcnow().hour
        mins = datetime.datetime.utcnow().minute
        second = datetime.datetime.utcnow().second
        self.calculate(y, m, d, h, mins, second)
        #self.print()

    def printnow(self):
        self.now()
        self.print()

    def listnow(self):
        """
        Returns current calculated values for the moon position
        as an array in the following format:
        [day, right_ascension, declination, latitude_info, longitude_info, azimuth,
         altitude, phase_info, mean_right_ascension, topocentric_declination, mpar?,
         msd, gmsto]
        """
        self.now()
        return self.list()

    def riseset(self, dis, utcdis):
        moonarray = self.risesetlist(dis, utcdis)
        print("moonrise time is {}:{}:{}".format(moonarray[0], moonarray[1], moonarray[2]))
        print("moonset time is {}:{}:{}".format(moonarray[3], moonarray[4], moonarray[5]))
        print("{} iterations".format(moonarray[6]))
        print("moon utc is {}".format(moonarray[7]))
        print("ra is {}".format(degrees(self.mra)))
        """
        y = datetime.datetime.utcnow().year
        m = datetime.datetime.utcnow().month
        d = datetime.datetime.utcnow().day + dis
        h = datetime.datetime.utcnow().hour
        mins = datetime.datetime.utcnow().minute
        second = datetime.datetime.utcnow().second
        self.calculate(y, m, d, h, mins, second)
        hmm = radians(-0.583) - self.msd
        mlha = (sin(hmm) - sin(radians(self.localLat)) * sin(self.mdecl))/(cos(radians(self.localLat)) * cos(self.mdecl))
        mlha = acos(mlha) * 12 * 15.0 / (15.04107 * pi)
        utcmoon = (degrees(self.mra) - self.gmsto * 15 - self.localLong) / 15
        a = mlha
        c = 1
        ff = 0
        while abs(c) > 0.000001:
            y = datetime.datetime.utcnow().year
            m = datetime.datetime.utcnow().month
            d = datetime.datetime.utcnow().day + dis
            h = utcmoon + a
            self.calculate(y, m, d, h, 0, 0)
            hmm = radians(-0.583) - self.msd
            mlha = (sin(hmm) - sin(radians(self.localLat)) * sin(self.mdecl))/(cos(radians(self.localLat)) * cos(self.mdecl))
            mlha = acos(mlha) * 12 * 15.0 / (15.04107 * pi)
            utcmoon = (degrees(self.mra) - self.gmsto * 15 - self.localLong) / 15
            b = mlha
            c = abs(b - a)
            a = mlha
            ff += 1
            if ff >= 100:
                break
            
        moonset = mlha + utcmoon
        c = 1
        
        while abs(c) > 0.000001:
            y = datetime.datetime.utcnow().year
            m = datetime.datetime.utcnow().month
            d = datetime.datetime.utcnow().day + dis
            h = utcmoon - a
            self.calculate(y, m, d, h, 0, 0)
            hmm = radians(-0.583) - self.msd
            mlha = (sin(hmm) - sin(radians(self.localLat)) * sin(self.mdecl))/(cos(radians(self.localLat)) * cos(self.mdecl))
            mlha = acos(mlha) * 12 * 15.0 / (15.04107 * pi)
            utcmoon = (degrees(self.mra) - self.gmsto * 15 - self.localLong) / 15
            b = mlha
            c = b - a
            a = mlha
            ff += 1
            if ff >= 200:
                break

        moonrise = utcmoon - mlha    
        moonrise = timeadjust(moonrise, utcdis)
        moonrise2 = decimaltominsecs(moonrise)
        moonset = timeadjust(moonset, utcdis)
        moonset2 = decimaltominsecs(moonset)"""
        

    def risesetlist(self, dis, utcdis):
        """
        Given a day displacement (dis) and utcdisplacement (utcdis),
        This function calculates and returns the corresponding moonset
        and moonrise. The result is given as an array with the following
        format:
        [moonrise_hour, moonrise_minute, moonrise_second, moonset_hour,
         moonset_minute, moonset_second, iterations_needed_calculate,
         moon_utc]
        """
        y = datetime.datetime.utcnow().year
        m = datetime.datetime.utcnow().month
        d = datetime.datetime.utcnow().day + dis
        h = datetime.datetime.utcnow().hour
        mins = datetime.datetime.utcnow().minute
        second = datetime.datetime.utcnow().second
        self.calculate(y, m, d, h, mins, second)
        hmm = radians(-0.583) - self.msd
        mlha = (sin(hmm) - sin(radians(self.localLat)) * sin(self.mdecl))/(cos(radians(self.localLat)) * cos(self.mdecl))
        mlha = acos(mlha) * 12 * 15.0 / (15.04107 * pi)
        utcmoon = (degrees(self.mra) - self.gmsto * 15 - self.localLong) / 15
        a = mlha
        c = 1
        ff = 0
        while abs(c) > 0.000001:
            y = datetime.datetime.utcnow().year
            m = datetime.datetime.utcnow().month
            d = datetime.datetime.utcnow().day + dis
            h = utcmoon + a
            self.calculate(y, m, d, h, 0, 0)
            hmm = radians(-0.583) - self.msd
            mlha = (sin(hmm) - sin(radians(self.localLat)) * sin(self.mdecl))/(cos(radians(self.localLat)) * cos(self.mdecl))
            mlha = acos(mlha) * 12 * 15.0 / (15.04107 * pi)
            utcmoon = (degrees(self.mra) - self.gmsto * 15 - self.localLong) / 15
            b = mlha
            c = abs(b - a)
            a = mlha
            ff += 1
            if ff >= 100:
                break
            
        moonset = mlha + utcmoon
        c = 1
        
        while abs(c) > 0.000001:
            y = datetime.datetime.utcnow().year
            m = datetime.datetime.utcnow().month
            d = datetime.datetime.utcnow().day + dis
            h = utcmoon - a
            self.calculate(y, m, d, h, 0, 0)
            hmm = radians(-0.583) - self.msd
            mlha = (sin(hmm) - sin(radians(self.localLat)) * sin(self.mdecl))/(cos(radians(self.localLat)) * cos(self.mdecl))
            mlha = acos(mlha) * 12 * 15.0 / (15.04107 * pi)
            utcmoon = (degrees(self.mra) - self.gmsto * 15 - self.localLong) / 15
            b = mlha
            c = b - a
            a = mlha
            ff += 1
            if ff >= 200:
                break

        moonrise = utcmoon - mlha    
        moonrise = timeadjust(moonrise, utcdis)
        moonrise2 = decimaltominsecs(moonrise)
        moonset = timeadjust(moonset, utcdis)
        moonset2 = decimaltominsecs(moonset)

        return [floor(moonrise), moonrise2[0], moonrise2[1], floor(moonset), moonset2[0], moonset2[1], ff, utcmoon]
        
    def print(self):
        #self.now()
        print("days: {}".format(self.d))
        print(self.ra)
        print("RA: " + str(degrees(self.ra)/15))
        print("Dec " + str(degrees(self.dec)))
        print(self.lstring)
        print(self.lonstring)
        print("Azimuth is {} degrees".format(str(degrees(self.az))))
        print("Corrected Alt is {} degrees".format(str(degrees(self.alttopoc))))
        print(self.phasestr)


