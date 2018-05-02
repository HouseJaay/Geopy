import obspy
from glob import glob
from collections import namedtuple
import matplotlib.pyplot as plt
import cartopy.crs as ccrs


class Event(object):
    def __init__(self, tr, filename):
        Stats = namedtuple('Stats', ['time', 'lat', 'lon', 'depth'])
        time = tr.stats.starttime + tr.stats.sac['o']
        self.evtinfo = Stats(time=time, lat=tr.stats.sac['evla'],
                             lon=tr.stats.sac['evlo'], depth=tr.stats.sac['evdp'])
        self.sta = {}
        if tr.stats.station in self.sta:
            raise NameError("duplicate event in %s" % tr.stats.station)
        else:
            self.sta[tr.stats.station] = filename

    def addstation(self, filename):
        station = (filename.split('/')[-1]).split('.')[1]
        self.sta[station] = filename

    def getfile(self, staname):
        return self.sta[staname]

    def getlatlon(self):
        return self.evtinfo.lat, self.evtinfo.lon

    def gettime(self):
        return self.evtinfo.time

    def __str__(self):
        temp = "%.2f  %.2f  %.2f  %d  " % (
            self.evtinfo.lat, self.evtinfo.lon, self.evtinfo.depth, len(self.sta))
        return temp + str(self.evtinfo.time)
        

class Events(object):
    def __init__(self):
        self.evts = {}

    def addfromdir(self, directory, wild):
        times = map(lambda x: x.split('/')[-1], glob(directory + wild))
        for time in times:
            files = glob(directory + time + '/*.Z')
            filename = files[0]
            tr = obspy.read(filename)[0]
            self.evts[time] = Event(tr, filename)
            for file in files:
                self.evts[time].addstation(file)

    def __iter__(self):
        return iter(self.evts.values())

    def __len__(self):
        return len(self.evts)

    def plotlocation(self):
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(
            1, 1, 1, projection=ccrs.AzimuthalEquidistant(105, 33))

        ax.set_global()
        ax.stock_img()
        ax.coastlines()

        for evt in self.evts.values():
            ax.plot(
                *evt.getlatlon()[::-1], marker='*', color='red',
                markersize=5, transform=ccrs.Geodetic())

        plt.show()


if __name__ == '__main__':
    directory = '/home/haosj/data/neTibet/data/'
    evts = Events()
    evts.addfromdir(directory, '2013*')
    for evt in evts:
        print(evt)
