import obspy
from glob import glob
from collections import namedtuple
from distaz import distaz
from scipy import signal
import numpy as np
import os


class Event(object):
    def __init__(self, tr, filename, name):
        self.PRANGE = (10, 80)
        self.pairs = {}
        self.name = name
        Stats = namedtuple('Stats', ['time', 'lat', 'lon', 'depth', 'mag'])
        time = tr.stats.starttime + tr.stats.sac['o']
        self.evtinfo = Stats(time=time, lat=tr.stats.sac['evla'],
                             lon=tr.stats.sac['evlo'], depth=tr.stats.sac['evdp'],
                             mag=tr.stats.sac['mag'])
        self.sta = {}
        if tr.stats.station in self.sta:
            raise NameError("duplicate event in %s" % tr.stats.station)
        else:
            self.sta[tr.stats.station] = filename

    def add_pair(self, pairname, dist):
        self.pairs[pairname] = dist

    def do(self, outdir, delta=1.0):
        if not self.pairs:
            return
        filtered = {}
        epic_center_dist = {}
        for s in self.sta:
            epic_center_dist[s] = obspy.read(
                self.sta[s], headonly=True)[0].stats.sac['dist']

        def def_b(x):
            return signal.firwin(
                1001, [1.0/(x+0.2), 1.0/(x-0.2)],
                window=('kaiser', 9), nyq=1/delta/2, pass_zero=False)
        B = list(map(def_b, range(*self.PRANGE)))
        epic_dist_max = max(epic_center_dist.values())
        epic_dist_min = min(epic_center_dist.values())
        starttime = self.gettime() + epic_dist_min/6.0
        endtime = self.gettime() + epic_dist_max/2.0
        npts = int((endtime - starttime) / delta)
        for s in self.sta:
            # TODO optimize space consume
            st = obspy.read(self.sta[s])
            if all(np.isnan(st[0].data)):
                continue
            try:
                data_cut = cut(st, starttime, npts)
            except IndexError:
                print(self.name+" "+s+" cut error")
                continue
            filtered[s] = norm2d(multi_filter(B, data_cut))
        for pair in self.pairs:
            outpath = outdir + pair + '/'
            if not os.path.exists(outpath):
                os.mkdir(outpath)
            sta1, sta2 = pair.split('_')
            if not(sta1 in filtered and sta2 in filtered):
                continue
            if epic_center_dist[sta2] > epic_center_dist[sta1]:
                sta1, sta2 = sta2, sta1
            result = two_station(
                filtered[sta1], filtered[sta2], self.pairs[pair], delta, (2.5, 4.5))
            np.savetxt(outpath+self.name, result)

    def addstation(self, filename):
        station = (filename.split('/')[-1]).split('.')[1]
        self.sta[station] = filename

    def getfile(self, staname):
        return self.sta[staname]

    def getlatlon(self):
        return self.evtinfo.lat, self.evtinfo.lon

    def gettime(self):
        return self.evtinfo.time

    def getmag(self):
        return self.evtinfo.mag

    def getdepth(self):
        return self.evtinfo.depth

    def __str__(self):
        temp = "%.2f  %.2f  %.2f  %d  " % (
            self.evtinfo.lat, self.evtinfo.lon, self.evtinfo.depth, len(self.sta))
        return temp + str(self.evtinfo.time)


class Events(object):
    def __init__(self):
        self.evts = {}

    def addfromdir(self, directory, timewild, filewild='*.Z'):
        times = map(lambda x: x.split('/')[-1], glob(directory + timewild))
        for time in times:
            files = glob(directory + time + '/' + filewild)
            if len(files) > 0:
                filename = files[0]
                tr = obspy.read(filename)[0]
                self.evts[time] = Event(tr, filename, time)
                for file in files:
                    self.evts[time].addstation(file)

    def prep_pairs(self, pairs):
        for pair in pairs:
            pair.setevents(self)

    def __iter__(self):
        return iter(self.evts.values())

    def __len__(self):
        return len(self.evts)


class Pairs(object):
    def __init__(self, stafile):
        self.pairs = {}
        statemp = []
        with open(stafile, 'r') as f:
            for line in f:
                l = line[:-1].split(' ')
                statemp.append([l[0], float(l[1]), float(l[2])])
        for i in range(len(statemp)-1):
            for j in range(i+1, len(statemp)):
                key = self.getkey(statemp[i][0], statemp[j][0])
                self.pairs[key] = Pair(*statemp[i], *statemp[j])

    @staticmethod
    def getkey(name1, name2):
        if name1 < name2:
            key = name1 + '_' + name2
        else:
            key = name2 + '_' + name1
        return key

    def __iter__(self):
        return iter(self.pairs.values())


class Pair(object):
    """
    usage:
        pair.setevents -> pair.do_ts
        or pair.setevents -> pair.load
    """
    def __init__(self, sta1, lat1, lon1, sta2, lat2, lon2):
        self.name = Pairs.getkey(sta1, sta2)
        self.sta1 = sta1
        self.sta2 = sta2
        self.latlon1 = (lat1, lon1)
        self.latlon2 = (lat2, lon2)
        self.staAz = distaz(*self.latlon2, *self.latlon1).getAz()
        self.staDist = distaz(*self.latlon2, *self.latlon1).degreesToKilometers()
        self.evts = []
        self.dispfile = []
        self.disp = []
        self.PRANGE = (10, 80)

    def getdiffaz(self, event):
        return abs(self.staAz - distaz(*event.getlatlon(), *self.latlon1).getAz())

    def setevents(self, events):
        for evt in events:
            diff = self.getdiffaz(evt)
            if evt.getdepth() > 50:
                continue
            if distaz(*evt.getlatlon(), 36, 105).getDelta() > 120:
                continue
            if not(diff < 2.0 or abs(diff - 180) < 2.0):
                continue
            if evt.getmag() <= 5.5:
                continue
            dist = min(distaz(*evt.getlatlon(), *self.latlon1).getDelta(),
                       distaz(*evt.getlatlon(), *self.latlon2).getDelta())
            if dist < 10:
                continue
            try:
                evt.getfile(self.sta1)
                evt.getfile(self.sta2)
            except KeyError:
                continue
            self.evts.append(evt)
            evt.add_pair(self.name, self.staDist)
        return


def cut(st, start, npts):
    delta = st[0].stats.delta
    n0 = int((start - st[0].stats.starttime) / delta)
    if n0 < 0 or n0+npts >= len(st[0].data):
        raise IndexError
    return st[0].data[n0:n0+npts]


def multi_filter(B, data):
    mat = []
    for b in B:
        mat.append(signal.lfilter(b, 1, data))
    return np.array(mat)


def norm2d(mat):
    ma, mi = mat.max(axis=1), mat.min(axis=1)
    m = np.c_[abs(ma), abs(mi)].max(axis=1)
    return mat / m.reshape(len(m), 1)


def pick(cor, uini, u):
    j = 0
    for i in range(len(u)):
        if u[i] <= uini:
            j = i
            break
    if j == 0 or j == (len(u)-1):
        return -1
    if cor[j+1] > cor[j]:
        while j < (len(u)-1) and cor[j+1] > cor[j]:
            j += 1
        i = j
    elif cor[j-1] > cor[j]:
        while j > 0 and cor[j-1] > cor[j]:
            j -= 1
        i = j
    return u[i]


def two_station(data1, data2, dist, delta, vrange):
    npts = len(data1[0])
    len_cor = 2*npts - 1
    t = np.arange(1, int((len_cor + 1) / 2)) * delta
    v = dist / t
    mask = (v > vrange[0]) * (v < vrange[1])
    v = v[mask]
    COR = []
    for i in range(len(data1)):
        COR.append(signal.correlate(data1[i], data2[i], mode='full'))
    COR = np.array(COR)
    COR = COR[:, int((len_cor+1)/2):len_cor]
    COR = COR[:, mask]
    result = np.empty(len(data1))
    result[:] = np.nan
    index = 50
    uini = pick(COR[index], 4.0, v)
    result[index] = uini
    utemp = uini
    for i in range(index+1, len(COR)):
        utemp = pick(COR[i], utemp, v)
        if utemp > 0:
            result[i] = utemp
        else:
            break
    utemp = uini
    for i in range(index-1, -1, -1):
        utemp = pick(COR[i], utemp, v)
        if utemp > 0:
            result[i] = utemp
        else:
            break
    return result


def test_ts(f1, f2, refdisp):
    # dist f1 > dist f2
    import matplotlib.pyplot as plt
    dire = '/home/haosj/data/fk_data2/'
    # dire = '/home/haosj/work/Geopy/testdata/twostation/'
    st = obspy.read(dire+f1)
    st += obspy.read(dire+f2)
    delta = st[0].stats.delta

    def cut2(st, start, end):
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

    cut2(st, max(st[0].stats.starttime, st[1].stats.starttime),
         min(st[0].stats.endtime, st[1].stats.endtime))

    def def_b(x):
        return signal.firwin(
            1001, [1.0 / (x + 0.2), 1.0 / (x - 0.2)],
            window=('kaiser', 9), nyq=1 / delta / 2, pass_zero=False)
    B = list(map(def_b, range(10, 80)))
    mat1 = norm2d(multi_filter(B, st[0].data))
    mat2 = norm2d(multi_filter(B, st[1].data))
    result = two_station(
        mat1, mat2, st[0].stats.sac['dist']-st[1].stats.sac['dist'], delta, (2, 6))
    forward = np.loadtxt(dire + refdisp).transpose()
    fig, ax = plt.subplots()
    ax.set_xlim(10, 80)
    ax.set_ylim(3, 4.5)
    ax.plot(np.arange(10, 80), result, '.', color='blue')
    ax.plot(forward[0], forward[1], color='orange')
    plt.show()
    return result


if __name__ == '__main__':
    dire = '/home/haosj/data/neTibet/'
    pairs = Pairs(dire + 'sta_II.lst')
    evts = Events()
    evts.addfromdir(dire+'data/', '2015*')
    evts.prep_pairs(pairs)
    evt = evts.evts['20150102.002.082155.900']
    # result = test_ts('g70.z', 'g60.z', 'forward_rayleigh')
