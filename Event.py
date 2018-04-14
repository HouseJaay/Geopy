import obspy
from glob import glob
from collections import namedtuple

class Event(object):
    def __init__(self,tr,filename):
        Stats = namedtuple('Stats',['time','lat','lon','depth'])
        time = tr.stats.starttime + tr.stats.sac['o']
        self.evtinfo = Stats(time=time,lat=tr.stats.sac['evla'],
        lon=tr.stats.sac['evlo'],depth=tr.stats.sac['evdp'])
        self.sta = {}
        if tr.stats.station in self.sta:
            raise NameError("duplicate event in %s" % tr.stats.station)
        else:
            self.sta[tr.stats.station] = filename

    def addstation(self,tr,filename):
        if tr.stats.station in self.sta:
            raise NameError("duplicate event in %s" % tr.stats.station)
        else:
            self.sta[tr.stats.station] = filename

    def getfile(self,staname):
        return self.sta[staname]

    def getlatlon(self):
        return (self.evtinfo.lat, self.evtinfo.lon)

    def gettime(self):
        return self.evtinfo.time

    def __str__(self):
        temp = "%.2f  %.2f  %.2f  %d  " % (self.evtinfo.lat, 
        self.evtinfo.lon, self.evtinfo.depth, len(self.sta))
        return temp + str(self.evtinfo.time)
        

class Events(object):
    def __init__(self):
        self.evts = {}

    def addevents(self,files):
        for filename in files:
            tr = obspy.read(filename, headonly=True)[0]
            time = tr.stats.starttime + tr.stats.sac['o']
            name = "%.2f%.2f%.2f%d%d" % (
                tr.stats.sac['evla'],
                tr.stats.sac['evlo'], tr.stats.sac['evdp'],
                time.minute,time.second)
            if name in self.evts:
                self.evts[name].addstation(tr, filename)
            else:
                self.evts[name] = Event(tr, filename)

    def __iter__(self):
        return iter(self.evts.values())



if __name__ == '__main__':
    #dirs = glob('/home/haosj/data/tibet/ped/*')
    dirs = ['/home/haosj/data/tibet/ped/15639',
    '/home/haosj/data/tibet/ped/15640']
    evts = Events()
    for dire in dirs:
        files = glob(dire+'/X2.*.Z')
        evts.addevents(files)
    for evt in evts:
        print(evt)
