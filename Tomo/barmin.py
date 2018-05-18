import os
import numpy as np
from distaz import distaz
import math


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


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    data = '/home/haosj/data/neTibet/'
    ray, time, mis = read_dispersion(data+'result2/', data+'metadata/sta_36_south.lst', 10)
    rnodes = [(x, y) for x in range(30, 35) for y in range(100, 105)]
    F = define_F([1, 2], [2, 3], rnodes)
    plt.matshow(F)
    plt.colorbar()
    plt.show()
