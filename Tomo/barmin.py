import os
import numpy as np
from distaz import distaz
import math
import matplotlib.pyplot as plt
import cartopy.crs as ccrs


def read_dispersion(disp_dir, stationlst, index):
    """
    read dispersion file, return start and end position of ray path,
        travel time and corresponding misfit
    :param disp_dir: directory of dispersion file
    :param stationlst: station info file
    :param index: correspond to certain peroid
    :return: (ray, time, mis)
        ray: lat1, lon1, lat2, lon2
        time: travel time
        mis: misfit
    """
    ray, time, mis = [], [], []
    sta = {}
    with open(stationlst, 'r') as f:
        for line in f:
            l = line[:-1].split(' ')
            sta[l[0]] = (float(l[1]), float(l[2]))
    disp_files = os.listdir(disp_dir)
    for disp_file in disp_files:
        disp = np.loadtxt(disp_dir+disp_file)
        if sum(np.isnan(disp[:, index])) == 0:
            vel, mis_v = disp[0, index], disp[1, index]
            sta1, sta2 = disp_file.split('_')
            latlon1, latlon2 = sta[sta1], sta[sta2]
            distance = distaz(*latlon2, *latlon1).degreesToKilometers()
            time.append(vel*distance)
            mis.append(mis_v*distance)
            ray.append(latlon1+latlon2)
    return ray, time, mis


def define_F(alpha, delta, rnodes):
    """
    define smooth matrix F
    :param alpha: smooth weight
    :param delta: smooth spatial width,
        len(delta)==len(alpha)==n+1
    :param rnodes: position of nodes
    :return: F matrix, (n+1)M x (n+1)M
    """
    def S(r1, r2, delta):
        return math.exp(-((r1[0]-r2[0])**2 + (r1[1]-r2[1])**2) / (2*delta**2))
    n, M = len(alpha)-1, len(rnodes)
    F = np.zeros([(n+1)*M, (n+1)*M], dtype='float')
    for k in range(len(alpha)):
        for j1 in range(k*M, (k+1)*M):
            for j2 in range(j1+1, (k+1)*M):
                s = S(rnodes[j1-k*M], rnodes[j2-k*M], delta[k])
                F[j1, j2] = s
                F[j2, j1] = s
            F[j1] /= np.sum(F[j1])
        for j in range(k*M, (k+1)*M):
            F[j, j] = -1
        F[k*M:(k+1)*M, k*M:(k+1)*M] *= -alpha[k]
    return F


def sphere2cartesian(colat, lon):
    theta = colat * math.pi / 180.0
    phi = lon * math.pi / 180.0
    z = math.cos(theta)
    x = math.sin(theta) * math.cos(phi)
    y = math.sin(theta) * math.sin(phi)
    return np.array([x, y, z])


def ray_dist_to_node(ray, rnode):
    """
    compute distance from point to great circle
    :param ray: start,end location, (colat1,lon1), (colat2,lon2)
    :param rnode: location of node
    :return: (bool, dist)
        if node is in range of arc, True
        distance in degrees
    """
    ray = [(90-ray[0][0], 180-ray[0][1]), (90-ray[1][0], 180-ray[1][1])]
    rnode = [90-rnode[0], 180-rnode[1]]
    a, b = sphere2cartesian(*ray[0]), sphere2cartesian(*ray[1])
    n = np.cross(a, b) / np.linalg.norm(np.cross(a, b))
    r = sphere2cartesian(*rnode)
    distance = abs(math.asin(np.dot(n, r))) * 180/math.pi
    rplane = r - np.dot(r, n) * n
    rplane = rplane / np.linalg.norm(rplane)
    in_arc = (
                abs(math.acos(np.dot(a, b)) -
                    (math.acos(np.dot(a, rplane)) + math.acos(np.dot(b, rplane))))
                < 0.01)
    return in_arc, distance


def test_ray_dist_to_node():
    from random import random
    from distaz import distaz
    ray = [
            ((random()-0.5)*180, (random()-0.5)*360), ((random()-0.5)*180, (random()-0.5)*360)]
    ax = plt.axes(projection=ccrs.Robinson())
    ax.set_global()
    ax.coastlines()
    plt.plot([ray[0][1], ray[1][1]], [ray[0][0], ray[1][0]], color='red', transform=ccrs.Geodetic())
    for lon in range(-180, 180, 10):
        for lat in range(-90, 90, 10):
            color = 'grey'
            inarc, dist = ray_dist_to_node(ray, (lat, lon))

            if dist < 5 and inarc:
                color = 'red'
            elif dist < 5:
                d1 = distaz(lat, lon, *ray[1]).getDelta()
                d2 = distaz(lat, lon, *ray[0]).getDelta()
                if d1 < 5 or d2 < 5:
                    color = 'red'

            plt.plot([lon], [lat], color=color, marker='.', transform=ccrs.Geodetic())
    plt.show()


def compute_ray_density(rays, rnodes, radius):
    """
    compute ray density to determine H
    :param rays: location of station pairs
    :param rnodes: location of nodes
    :param radius: radius of node
    :return: rho and chi, ray density and azimuthal distribution
    """
    rho = [0 for _ in range(len(rnodes))]
    azs = [[] for _ in range(len(rnodes))]
    for ray in rays:
        az = distaz(*ray[0], *ray[1]).getAz()
        if az > 180:
            az -= 180
        for i in range(len(rnodes)):
            inarc, dist = ray_dist_to_node(ray, rnodes[i])
            if dist < radius:
                d1 = distaz(*rnodes[i], *ray[0]).getDelta()
                d2 = distaz(*rnodes[i], *ray[1]).getDelta()
                if inarc or (d1 < radius or d2 < radius):
                    rho[i] += 1
                    azs[i].append(az)


if __name__ == '__main__':
    data = '/home/haosj/data/neTibet/'
    ray, time, mis = read_dispersion(data+'result2/', data+'metadata/sta_36_south.lst', 10)
    rnodes = [(x, y) for x in range(30, 35) for y in range(100, 105)]
    F = define_F([1, 2], [2, 3], rnodes)
    plt.matshow(F)
    plt.colorbar()
    plt.show(block=True)
    test_ray_dist_to_node()
