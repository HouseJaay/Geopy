import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from Geopy.Tomo.barmin import *


def number_of_path(indir, prange):
    files = glob(indir + '*')
    result = []
    for f in files:
        result.append(np.loadtxt(f)[0])
    result = np.array(result)
    num = np.sum(~np.isnan(result), axis=0)
    fig, ax = plt.subplots()
    ax.set_xlabel('Period(s)')
    ax.set_ylabel('Number of paths')
    ax.plot(np.arange(*prange), num)
    plt.show()


def read_ray_from_into(intofile):
    """
    read ray from intofile of DS2004
    :param intofile: intofile path
    :return: ray (lat1, lon1, lat2, lon2)
    """
    f = open(intofile, 'r')
    lines = f.readlines()
    rays = []
    for i in range(1, len(lines), 4):
        rays.append(list(map(lambda x: float(x), lines[i].split())))
    return rays


def get_density(lat1, lat2, lon1, lon2, radius, intofile, threshold):
    """
    use ray density filter tomography
    :return: dict key: lat_lon value:boolean True means keep corresponding node
    """
    rays = read_ray_from_into(intofile)
    rnodes = [(x, y) for x in np.arange(lat1, lat2+radius, radius)
              for y in np.arange(lon1, lon2+radius, radius)]
    rho, chi = compute_ray_density(rays, rnodes, radius)
    density = set()
    for i in range(len(rnodes)):
        if rho[i] > threshold:
            density.add(str(rnodes[i][0]) + '_' + str(rnodes[i][1]))
    return density


if __name__ == '__main__':
    number_of_path('/home/haosj/data/last2/result2_rejected2/', (10, 80))
    rays = read_ray_from_into('/home/haosj/data/last2/tomo/tibet_into')
