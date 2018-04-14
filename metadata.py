import distaz
import matplotlib.pyplot as plt
import numpy as np
import pandas
from obspy.core import UTCDateTime


# cut every seismogram in a stream
def cut(st, start, end):
    delta = st[0].stats.delta
    width = end - start
    npts = int(width/delta)
    # check
    if(start > end):
        raise ValueError('starttime must be smaller than endtime')
    for i in range(len(st)):
        if(st[i].stats.delta != delta):
            raise ValueError('all data must have same delta')
        if(st[i].stats.starttime > start or st[i].stats.endtime < end):
            raise ValueError('invalid time range')
    for i in range(len(st)):
        n = int((start - st[i].stats.starttime)/delta)
        st[i].data = st[i].data[n:n+npts]
        st[i].stats.starttime = start

# plot seismogram in a stream
def plot(st):
    n = st[0].stats.npts
    delta = st[0].stats.delta
    t = np.arange(n)*delta
    rows = len(st)
    fig,axes = plt.subplots(nrows=rows,ncols=1)
    for i in range(rows):
        axes[i].plot(t,st[i].data,color='black')
    plt.show()

# input file is bqmail output format
def read_station(filename):
    """
    input file format:
    net station lat lon start end
    """
    col_name = ['net','station','lat','lon','start','end']
    sta = pandas.read_table(filename,sep='\s+',names=col_name)
    sta['start'] = sta['start'].map(read_date)
    sta['end'] = sta['end'].map(read_date)
    return sta

# read date in year/month/day to UTCDateTime format
def read_date(date):
    """
    read date year/month/day
    return UTCDateTime
    """
    temp = date.split('/')
    year = int(temp[0])
    month = int(temp[1])
    day = int(temp[2])
    return UTCDateTime(year,month,day)

def merge_date(start1,end1,start2,end2):
    if(start1 > start2):
        start = start1
    else:
        start = start2
    if(end1 > end2):
        end = end2
    else:
        end = end1
    if(end < start):
        raise ValueError('station pair must have common work date')
    else:
        return start,end

# input value is output of read_table
def mk_sta_pairs(sta):
    col_name = ['net1','station1','lat1','lon1','net2','station2','lat2','lon2','start','end']
    pair = pandas.DataFrame(columns=col_name)
    row = 0
    for i in range(len(sta)-1):
        for j in range(i+1,len(sta)):
            temp1 = [sta['net'][i],sta['station'][i],sta['lat'][i],sta['lon'][i]]
            temp2 = [sta['net'][j],sta['station'][j],sta['lat'][j],sta['lon'][j]]
            try:
                start,end = merge_date(sta['start'][i],sta['end'][i],sta['start'][j],sta['end'][j])
            except ValueError:
                print("date error:%s %s and %s %s" % (sta['net'][i],sta['station'][i],sta['net'][j],sta['station'][j]))
                continue
            pair.loc[row] = temp1 + temp2 + [start,end]
            row += 1
    return pair

def cal_dist(row):
    d = distaz.distaz(row['lat1'],row['lon1'],row['lat2'],row['lon2']).getDelta()
    return d

# input file in bqmail output format
def read_event(filename):
    """
    input file format:
    year month day jday hour min sec lat lon dep mw
    """
    col_name = ['year','month','day','jday','hour','min','sec','lat','lon','dep','mw']
    evt = pandas.read_table(filename, sep='\s+', names=col_name)
    time = lambda x : UTCDateTime(int(x['year']),int(x['month']),int(x['day']),int(x['hour']),int(x['min']),int(x['sec']))
    evt['time'] = evt.apply(time,axis='columns')
    del evt['year'], evt['month'], evt['day'], evt['jday'], evt['hour'], evt['min'], evt['sec']
    return evt

# input list of filenames
def read_data(filelist):
    pass
def plot_data(data):
    pass
