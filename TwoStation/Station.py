from distaz import distaz
from Geopy.Event import Events,Event
from glob import glob
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from Geopy.TwoStation import two_station as ts
import obspy


class Pair(object):
    def __init__(self, sta1, lat1, lon1, sta2, lat2, lon2):
        self.sta1 = sta1
        self.sta2 = sta2
        self.latlon1 = (lat1, lon1)
        self.latlon2 = (lat2, lon2)
        self.staAz = distaz(*self.latlon2, *self.latlon1).getAz()
        self.staDist = distaz(*self.latlon2, *self.latlon1).degreesToKilometers()
        self.evts = []

    def getdiffaz(self, event):
        return abs(self.staAz - distaz(*event.getlatlon(), *self.latlon1).getAz())

    def setevents(self,events):
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
                *evt.getlatlon()[::-1],marker='*',color='red',
                markersize=5, transform=ccrs.Geodetic())

        plt.show()

    def _cut(self, st, start, end):
        delta = st[0].stats.delta
        width = end - start
        npts = int(width / delta)
        n0 = int((start - st[0].stats.starttime) / delta)
        n1 = int((start - st[1].stats.starttime) / delta)
        st[0].data = st[0].data[n0:n0+npts]
        st[1].data = st[1].data[n1:n1+npts]
        st[0].stats.starttime = start
        st[1].stats.starttime = start

    def do_ts(self):
        for evt in self.evts:
            st = obspy.read(evt.getfile(self.sta1))
            st += obspy.read(evt.getfile(self.sta2))
            dist = (st[0].stats.sac.dist + st[1].stats.sac.dist)/ 2.0
            starttime = evt.gettime() + dist/6
            endtime = evt.gettime() + dist/2
            self._cut(st, starttime, endtime)
            print(st)
            ts.two_station(st, self.staDist, (2, 6), (20, 80), '/home/haosj/data/outtemp')


if __name__ == '__main__':
    pair = Pair('15639', 38.681, 104.352, '61061', 36.526, 108.772)
    dirs = ['/home/haosj/data/tibet/ped/15639',
            '/home/haosj/data/tibet/ped/61061']
    evts = Events()
    for dire in dirs:
        files = glob(dire+'/X2.*.Z')
        evts.addevents(files)
    pair.setevents(evts)
    for evt in pair.evts:
        print(evt)
    pair.plotlocation()
    pair.do_ts()
