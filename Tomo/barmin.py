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
        ray = [(ray[0], ray[1]), (ray[2], ray[3])]
        az = distaz(*ray[0], *ray[1]).getAz()
        if az >= 180:
            az -= 180
        for i in range(len(rnodes)):
            inarc, dist = ray_dist_to_node(ray, rnodes[i])
            if dist < radius:
                d1 = distaz(*rnodes[i], *ray[0]).getDelta()
                d2 = distaz(*rnodes[i], *ray[1]).getDelta()
                if inarc or (d1 < radius or d2 < radius):
                    rho[i] += 1
                    azs[i].append(az)
                    # TODO fix bug. azimuths is not constant on one ray
    chi = [0 for _ in range(len(rnodes))]
    for i in range(len(chi)):
        if len(azs[i]) == 0:
            chi[i] = 0
        else:
            hist = [0 for _ in range(10)]
            for az in azs[i]:
                hist[int(az)//18] += 1
            chi[i] = sum(hist) / (10 * max(hist))
    return rho, chi


def define_H(rho, chi, beta, lambda0, thres_chi):
    """
    define regularization matrix H
    :param rho: ray path density
    :param chi: ray azimuth density
    :param beta: damp parameter
    :param lambda0: parameter define function from rho to damp coefficient
    :param thres_chi: if chi<thres_chi, don't inv anisotropy model
    :return: H  M*(n+1) x M*(n+1)
    """
    M, n = len(rho), len(beta) - 1
    H = np.zeros([M*(n+1), M*(n+1)])
    for i in range(M):
        H[i, i] = beta[0] * math.exp(-lambda0 * rho[i])
    for k in range(1, len(beta)):
        for j in range(M):
            H[k*M+j, k*M+j] = beta[k] * (0 if chi[j] < thres_chi else 1)
    return H


def plot_ray(rays, plotlim):
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent(plotlim, crs=ccrs.PlateCarree())
    ax.stock_img()
    ax.coastlines()
    for ray in rays:
        ax.plot([ray[1], ray[3]], [ray[0], ray[2]], color='grey', transform=ccrs.Geodetic())
    plt.show()


def plot_density(density, rnodes, r):
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent(r, crs=ccrs.PlateCarree())
    ax.stock_img()
    ax.coastlines()
    rnodes = np.array(rnodes)
    lats = sorted(list(set(rnodes[:, 0])))
    lons = sorted(list(set(rnodes[:, 1])))
    density = np.array(density).reshape(len(lats), len(lons))
    temp = ax.contourf(lons, lats, density, transform=ccrs.PlateCarree())
    fig.colorbar(temp)
    plt.show()


def gmt_density(val, rnodes, filename):
    f = open(filename, 'w')
    for i, rnode in enumerate(rnodes):
        f.write("%f %f %f\n" % (rnode[1], rnode[0], val[i]))
    f.close()


def gmt_ray(rays, filename):
    f = open(filename, 'w')
    for ray in rays:
        f.write(">\n%f %f\n%f %f\n" % (ray[1], ray[0], ray[3], ray[2]))
    f.close()


# project on a single plane, because current research area is small
# change later
def proj_to_cube(lat, lon, clon=105):
    """
    project sphere point on square, edge length = 2
    :param lat: latitude (-45 ~ 45)
    :param lon: longitude (clon-45 ~ clon+45)
    :param clon: central longitude
    :return: coordinate on plane
    """
    lat = lat * math.pi / 180.0
    lon = lon * math.pi / 180.0
    clon = clon * math.pi / 180.0
    y = math.tan(lat)
    x = math.tan(lon-clon)
    return x, y


def proj_to_sphere(x, y, clon=105):
    lat = math.atan(y) * 180 / math.pi
    lon = math.atan(x) * 180 / math.pi + clon
    return lat, lon


if __name__ == '__main__':
    data = '/home/haosj/data/neTibet/'
    laran, loran = (30, 37), (101, 109)
    radius = 0.5
    ray, time, mis = read_dispersion(
            data+'result3/', data+'metadata/sta_36_south.lst', 10)
    rnodes = [
            (x, y) for x in np.arange(laran[0], laran[1]+radius, radius)
            for y in np.arange(loran[0], loran[1]+radius, radius)]
    # F = define_F([1, 2], [2, 3], rnodes)
    # plt.matshow(F)
    # plt.colorbar()
    # plt.show(block=True)
    # test_ray_dist_to_node()
    # plot_ray(ray, loran+laran)
    # rho, chi = compute_ray_density(ray, rnodes, radius/2.0)
    # gmt_density(rho, rnodes, './plot/rho')
    # gmt_density(chi, rnodes, './plot/chi')
    # plot_density(rho, rnodes, loran+laran)
    gmt_ray(ray, './plot/ray')
