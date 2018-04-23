#! /usr/bin/env python

import obspy
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
import pandas
from obspy.core import UTCDateTime
import obspy.signal.filter

VRANGE = (3,5)
CH = 'BHZ'


def pick_global(array):
    arr_env = obspy.signal.filter.envelope(array)
    return arr_env.argmax()


def window(array,n0,width):
    l = n0 - width
    r = n0 + width
    if l<0: l=0
    if r>=len(array): r=len(array)-1
    w = np.zeros(len(array))
    w[l:r] = 1
    array = w*array
    return array


def norm(array):
    ma,mi = array.max(),array.min()
    m = max(abs(ma),abs(mi))
    return array/m


def pick(cor, uini, u):
    for i in range(len(u)):
        if u[i] <= uini:
            j = i
            break
    if j==0 or j==(len(u)-1):
        return -1
    if cor[j+1]>cor[j]:
        while j<(len(u)-1) and cor[j+1]>cor[j]:
            j+=1
        i=j
    elif cor[j-1]>cor[j]:
        while j>0 and cor[j-1]>cor[j]:
            j-=1
        i=j
    return u[i]    


def onclick(event):
    global click_x,click_y
    click_x,click_y=event.xdata,event.ydata


def two_station(st ,dist,vrange,prange):
    delta = st[0].stats.delta
    npts = st[0].stats.npts

    len_cor = 2*npts-1
    t = np.arange(1,int((len_cor+1)/2))*delta
    v = dist/t
    mask = (v>vrange[0]) * (v<vrange[1])
    v = v[mask]
    p = np.arange(prange[0], prange[1])
    COR = np.zeros((len(p),len(v)))
    V,P = np.meshgrid(v,p)
    
# sequence
    if(st[0].stats.sac.dist > st[1].stats.sac.dist):
        tr1 = st[0].copy()
        tr2 = st[1].copy()
    else:
        tr1 = st[1].copy()
        tr2 = st[0].copy()
# plot prepare
    rows = int((prange[1]-prange[0])/10)
    fig,axes = plt.subplots(nrows = rows+1, ncols=2) 
    ax_t = np.arange(len(tr1.data))*delta
    axes[0,0].plot(ax_t,tr1.data,'b-')
    axes[0,0].set_title(tr1.stats.station)
    axes[0,1].plot(ax_t,tr2.data,'b-')
    axes[0,1].set_title(tr2.stats.station)

    row=0
    for period in range(prange[0],prange[1]):
        b = signal.firwin(1001,[1.0/(period+0.2),1.0/(period-0.2)],window=('kaiser',9),nyq=1/delta/2,pass_zero=False)
# filter
        array1 = signal.lfilter(b,1,tr1.data)
        array2 = signal.lfilter(b,1,tr2.data)
# normalize
        array1 = norm(array1)
        array1 = signal.detrend(array1)
        array2 = norm(array2)
        array2 = signal.detrend(array2)
# window and plot
        #a1_temp = np.copy(array1)
        #a2_temp = np.copy(array2)
        #width = int(2*period/delta)
        #t1 = pick_global(array1)
        #t2 = pick_global(array2)
        #array1 = window(array1,t1,width)
        #array2 = window(array2,t2,width)
        if(row%10==0):
            nrow = int(row/10)+1
            #axes[nrow,0].plot(ax_t,a1_temp,'b-')
            #axes[nrow,1].plot(ax_t,a2_temp,'b-')
            axes[nrow,0].plot(ax_t,array1,'b-')
            axes[nrow,1].plot(ax_t,array2,'b-')

# correlate , first input signal has larger epicenter distance
        corr = signal.correlate(array1,array2,mode='full')
# data prepare
        cor = corr[int((len_cor+1)/2):len_cor]
        cor = cor[mask]
        cor = norm(cor)
        COR[row] = cor
        row += 1
# pick
    plt.show()

    fig,ax = plt.subplots()
    cf = ax.contourf(P,V,COR)
    fig.colorbar(cf)
    plt.show()
    result = np.empty(prange[1]-prange[0])
    result[:] = np.NaN
    pini = 60
    uini = pick(COR[pini-prange[0]], 4.0, v)

    result[pini-prange[0]] = uini
    utemp = uini
    for period in range(pini+1, prange[1], 1):
        utemp = pick(COR[period-prange[0]], utemp, v)
        if(utemp > 0):
            result[period-prange[0]] = utemp
        else:
            break
    utemp=uini
    for period in range(pini-1,prange[0]-1,-1):
        utemp = pick(COR[period-prange[0]], utemp, v)
        if(utemp > 0):
            result[period-prange[0]] = utemp
        else:
            break
    return result


def do_ts(Disp,PRANGE=(20,60)):
    for index,e in Disp.evt.iterrows():
        dist = abs(e['dist'][0]-e['dist'][1])
        two_station(e['data1'],e['data2'],dist,VRANGE,PRANGE,e['disp'])


if __name__ == '__main__':
    testdir = '/home/haosj/work/Geopy/testdata/twostation/'
    st = obspy.read(testdir + '*.z')
    print(st)
    result = two_station(st, 1111.95, (2, 6), (10, 100))
    forward = np.loadtxt(testdir + 'forward').transpose()
    fig, ax = plt.subplots()
    ax.set_xlim(10, 100)
    ax.set_ylim(3, 4)
    ax.plot(np.arange(10, 100), result)
    ax.plot(forward[0], forward[1])
    plt.show()