from pyproj import Geod
from Geopy.Tomo.barmin import sphere2cartesian
import numpy as np
from math import *
import random


def model1(lat, lon):
    """
    :param lat:
    :param lon:
    :return: [cI, A1, A2, A3, A4]
    """
    lat, lon = int(lat), int(lon)
    cI = ((lat % 2 and lon % 2) or (not lat % 2 and not lon % 2)) * 0.5 + 3.2
    return cI, 0, 0, 0, 0


def model_fake_anis(lat, lon):
    lat, lon = int(lat), int(lon)
    if lon < 105:
        cI = 3.2
    else:
        cI = 3.7
    return cI, 0, 0, 0, 0


def model2(lat, lon):
    """
    :param lat:
    :param lon:
    :return: [cI, A1, A2, A3, A4]
    """
    lat, lon = int(lat), int(lon)
    block = ((lat % 2 and lon % 2) or (not lat % 2 and not lon % 2))
    cI = block * 0.5 + 3.2
    la, lo = lat//2, lon//2
    block2 = ((la % 2 and lo % 2) or (not la % 2 and not lo % 2))
    A1 = block2 * 0.02 * cI
    A2 = (1 - block2) * 0.02 * cI
    return cI, A1, A2, 0, 0


def model3(lat, lon):
    """
    :param lat:
    :param lon:
    :return: [cI, A1, A2, A3, A4]
    """
    lat, lon = int(lat), int(lon)
    block = ((lat % 2 and lon % 2) or (not lat % 2 and not lon % 2))
    cI = block * 0.5 + 3.2
    la, lo = lat//2, lon//2
    block2 = ((la % 2 and lo % 2) or (not la % 2 and not lo % 2))
    A1 = 0.01 * cI * (-1)**(block2+1)
    A2 = 0
    return cI, A1, A2, 0, 0


def model_test(lat, lon):
    """
    :param lat:
    :param lon:
    :return: [cI, A1, A2, A3, A4]
    """
    lat, lon = int(lat), int(lon)
    block = ((lat % 2 and lon % 2) or (not lat % 2 and not lon % 2))
    cI = block * 0.5 + 3.2
    la, lo = lat//2, lon//2
    block2 = ((la % 2 and lo % 2) or (not la % 2 and not lo % 2))
    A2 = -0.01 * cI
    A1 = 0
    return cI, A1, A2, 0, 0


def write_model(model, ins, anis=None):
    f = open(ins, 'w')
    for lon in np.arange(-180, 180, 0.1):
        for lat in np.arange(-90, 90, 0.1):
            f.write("%f %f %f\n" % (lon, lat, model(lat, lon)[0]))
    f.close()
    if anis:
        f = open(anis, 'w')
        for lon in np.arange(-180, 180, 0.5):
            for lat in np.arange(-90, 90, 0.5):
                a1, a2 = model(lat, lon)[1:3]
                theta = atan2(a1, a2)
                theta = theta * 180 / pi
                azfast = (90 - theta) / 2
                gmt_azfast = 90 - azfast
                amp = sqrt(a1 ** 2 + a2 ** 2)
                gmt_amp = amp / (0.01*4)
                f.write("%f %f %f %fc 0.1c\n" % (lon, lat, gmt_azfast, gmt_amp))


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


def syn_ray(ray, model):
    """
    compute travel time along given ray
    :param ray: (lat1, lon1, lat2, lon2)
    :param model: function(lat, lon), return model
    :return: travel time, velocity
    """
    g = Geod(ellps='sphere')
    distance = g.inv(ray[1], ray[0], ray[3], ray[2])[2] / 1000.0
    points = g.npts(ray[1], ray[0], ray[3], ray[2], int(distance)//1)
    time = 0
    for i in range(len(points)-1):
        lo1, la1, lo2, la2 = points[i] + points[i+1]
        az, baz, d = g.inv(lo1, la1, lo2, la2)
        d /= 1000.0
        if az < 0:
            az += 180
        phi = az * pi / 180
        m = model(la1, lo1)
        c = m[0]
        c += m[1] * cos(2*phi) + m[2] * sin(2*phi) +\
            m[3] * cos(4*phi) + m[4] * sin(4*phi)
        time += d / c
    return time, distance / time


def checkboard_ds2004(into_file, model, out, gauss_std=0):
    """
    form intomodesVs type file for DS2004
    :param into_file: intomodesVs type file contain rays
    :param model: checkboard model
    :param out: output intomodesVs file
    :param gauss_std: add gauss noise
    :return:
    """
    fin = open(into_file, 'r')
    fout = open(out, 'w')
    lines = fin.readlines()
    for i in range(0, len(lines), 4):
        cur = lines[i:i+4]
        vel = syn_ray(list(map(lambda x: float(x), cur[1].split())), model)[1]
        vel += random.gauss(0, gauss_std)
        cur[2] = str(vel) + '\n'
        fout.write(''.join(cur))
    fin.close()
    fout.close()


if __name__ == '__main__':
    write_model(model_test, 'tomo_model_test', 'tomo_anvs_model_test')
    # test_integrate(0, 5, 80, 130)
    # t = syn_ray((30, 100, 35, 110), model1)
