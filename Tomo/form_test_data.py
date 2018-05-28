from pyproj import Geod
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from Geopy.Tomo.barmin import sphere2cartesian
import numpy as np


def model1(lat, lon):
    """
    :param lat:
    :param lon:
    :return: [cI, A1, A2, A3, A4]
    """
    cI = ((lat // 5) % 2 and (lon // 5) % 2) + 3
    return cI, 0, 0, 0, 0


def write_model(model, name):
    f = open(name, 'w')
    for lon in range(-180, 180):
        for lat in range(-90, 90):
            f.write("%f %f %f\n" % (lon, lat, model1(lat, lon)[0]))
    f.close()


def test_integrate(lat1, lon1, lat2, lon2):
    g = Geod(ellps='sphere')
    distance1 = g.inv(lon1, lat1, lon2, lat2)[2] / 1000.0
    points = g.npts(lon1, lat1, lon2, lat2, int(distance1)//1)
    distance2 = 0
    distance3 = 0
    res = 0
    for i in range(len(points)-1):
        lo1, la1, lo2, la2 = points[i] + points[i+1]
        p1 = sphere2cartesian(90-la1, 180-lo1) * 6371
        p2 = sphere2cartesian(90-la2, 180-lo2) * 6371
        distance2 += np.linalg.norm(p1-p2)
        az, baz, d = g.inv(lo1, la1, lo2, la2)
        distance3 += d
        res += abs(abs(az - baz) - 180)
    distance3 /= 1000.0
    res /= len(points)
    print(distance1, distance2, distance3)
    print(res)


if __name__ == '__main__':
    # write_model(model1, 'model1')
    test_integrate(0, 5, 80, 130)
