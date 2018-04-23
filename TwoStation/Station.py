from distaz import distaz
from Geopy.Event import Events
from glob import glob
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from Geopy.TwoStation import two_station as ts
import obspy
import numpy as np
from scipy import signal
import os


class Pair(object):
    def __init__(self, sta1, lat1, lon1, sta2, lat2, lon2):
        self.sta1 = sta1
        self.sta2 = sta2
        self.latlon1 = (lat1, lon1)
        self.latlon2 = (lat2, lon2)
        self.staAz = distaz(*self.latlon2, *self.latlon1).getAz()
        self.staDist = distaz(*self.latlon2, *self.latlon1).degreesToKilometers()
        self.evts = []
        self.dispfile = []
        self.disp = []
        self.PRANGE = (20, 80)

    def getdiffaz(self, event):
        return abs(self.staAz - distaz(*event.getlatlon(), *self.latlon1).getAz())

    def setevents(self, events):
        for evt in events:
            flag = True
            diff = self.getdiffaz(evt)
            if not(diff < 2.0 or abs(diff - 180) < 2.0):
                flag = False
            try:
                evt.getfile(self.sta1)
                evt.getfile(self.sta2)
            except KeyError:
                flag = False
            if flag:
                self.evts.append(evt)
        return

    def plotlocation(self):
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(
            1, 1, 1, projection=ccrs.AzimuthalEquidistant(*self.latlon1[::-1]))

        ax.set_global()
        ax.stock_img()
        ax.coastlines()
        
        ax.plot(*self.latlon1[::-1], marker='^', markersize=10, transform=ccrs.Geodetic())
        for evt in self.evts:
            ax.plot(
                *evt.getlatlon()[::-1], marker='*', color='red',
                markersize=5, transform=ccrs.Geodetic())

        plt.show()

    @staticmethod
    def cut(st, start, end):
        delta = st[0].stats.delta
        width = end - start
        npts = int(width / delta)
        n0 = int((start - st[0].stats.starttime) / delta)
        n1 = int((start - st[1].stats.starttime) / delta)
        if n0 < 0 or n1 < 0 or n0+npts >= len(st[0].data) or n1+npts >= len(st[1].data):
            raise IndexError
        st[0].data = st[0].data[n0:n0+npts]
        st[1].data = st[1].data[n1:n1+npts]
        st[0].stats.starttime = start
        st[1].stats.starttime = start

    def do_ts(self):
        snrthreshold = 7
        for evt in self.evts:
            st = obspy.read(evt.getfile(self.sta1))
            st += obspy.read(evt.getfile(self.sta2))
            dist = (st[0].stats.sac.dist + st[1].stats.sac.dist) / 2.0
            starttime = evt.gettime() + dist/6
            endtime = evt.gettime() + dist/2
            try:
                Pair.cut(st, starttime, endtime)
            except IndexError:
                continue
            if snr(st[0]) < snrthreshold or snr(st[1]) < snrthreshold:
                continue
            print(st)
            print(snr(st[0]), snr(st[1]))
            out_path = '/home/haosj/data/tibet/result/'
            out_path += self.sta1 + '_' + self.sta2 + '/'
            if not os.path.exists(out_path):
                os.mkdir(out_path)
            out_path += str(evt.gettime())
            phvel = ts.two_station(st, self.staDist, (2, 6), self.PRANGE)
            #TODO selection
            hand = input('input d to delete this result:')
            if hand != 'd':
                self.disp.append(phvel)
            #TODO write to file
            self.dispfile.append(out_path)
        self.disp = np.array(self.disp)

    def plot(self):
        fig, ax = plt.subplots()
        ax.set_xlabel('period(s)')
        ax.set_ylabel('phase velocity(km/s)')
        ax.set_xlim(*self.PRANGE)
        ax.set_ylim(2.9, 4.2)
        for i in range(len(self.disp)):
            plt.plot(np.arange(*self.PRANGE), self.disp[i])
        plt.show()


class Pairs(object):
    def __init__(self, stafile, events):
        self.events = events
        self.pairs = {}
        statemp = []
        with open(stafile, 'r') as f:
            for line in f:
                l = line[:-1].split(' ')
                statemp.append([l[0], float(l[1]), float(l[2])])
        for i in range(len(statemp)-1):
            for j in range(i+1, len(statemp)):
                key = self._getkey(statemp[i][0], statemp[j][0])
                self.pairs[key] = Pair(*statemp[i], *statemp[j])

    def getpair(self, sta1, sta2):
        key = self._getkey(sta1, sta2)
        pair = self.pairs[key]
        if pair.disp == []:
            pair.setevents(self.events)
            pair.plotlocation()
            pair.do_ts()
        return pair

    def _getkey(self, name1, name2):
        if name1 < name2:
            key = name1 + '_' + name2
        else:
            key = name2 + '_' + name1
        return key

    def plotdisp(self, *stanames):
        colors = ['b', 'r', 'g', 'c', 'm', 'y', 'k', 'w']
        PRANGE = (20, 80)
        if len(stanames) > len(colors):
            print('too many stations')
            return
        else:
            colors = iter(colors)
        fig, ax = plt.subplots()
        ax.set_xlabel('period(s)')
        ax.set_ylabel('phase velocity(km/s)')
        ax.set_xlim(*PRANGE)
        ax.set_ylim(2.9, 4.2)
        handle = []
        for staname in stanames:
            color = next(colors)
            key = self._getkey(*staname)
            pair = self.getpair(*staname)
            for i in range(len(pair.disp)):
                temp, = ax.plot(np.arange(*PRANGE),
                                pair.disp[i], '-', color=color, label=key)
                if i == 0:
                    handle.append(temp)
        ax.legend(handles=handle)
        plt.show()


def plotdispfile(filelist, pmin, pmax):
    fig, ax = plt.subplots()
    ax.set_xlabel('period(s)')
    ax.set_ylabel('phase velocity')
    ax.set_ylim(3, 5)
    ax.set_xlim(pmin, pmax)
    for f in filelist:
        disp = np.loadtxt(f)
        disp = disp.transpose()
        indexer = disp[0].argsort()
        ax.plot(disp[0][indexer], disp[1][indexer], '-')
    plt.show()


def snr(tr):
    data = abs(tr.data)
    mid = len(data)//2
    return data.max()/data[mid:].mean()


if __name__ == '__main__':
    dirs = glob('/home/haosj/data/tibet/ped/*')
    evts = Events()
    for dire in dirs:
        files = glob(dire+'/X2.*.Z')
        evts.addevents(files)
    # pair = Pair('15639', 38.681, 104.352, '61061', 36.526, 108.772)
    pairs = Pairs('/home/haosj/data/tibet/metadata/ordos_sta.lst', evts)
    pair = pairs.getpair('15639', '61061')