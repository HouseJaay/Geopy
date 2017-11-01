from scipy import signal
import obspy
import matplotlib.pyplot as plt 
import numpy as np
from Geopy import metadata as mt
from obspy.signal.filter import envelope

def cal_zhratio(zfile,rfile,freqs,bpwidth=0.002,outname=None,plot=None):
    try:
        st = obspy.read(zfile)
        st += obspy.read(rfile)
    except Exception:
        print("file error %s %s" %(zfile,rfile))
        return 0
    delta = st[0].stats.delta
    def_b = lambda f:signal.firwin(1001,[f-bpwidth,f+bpwidth],
    window=('kaiser',9),nyq=1/delta/2,pass_zero=False)
    B = list(map(def_b,freqs))
    Z = list(map(lambda b:signal.lfilter(b,1,st[0].data),B))
    H = list(map(lambda b:signal.lfilter(b,1,st[1].data),B))
    Z_env = list(map(envelope,Z))
    H_env = list(map(envelope,H))
    ratio = np.array(list(map(lambda e:max(e[0])/max(e[1]),
    zip(Z_env,H_env))))

    result = np.vstack((freqs,ratio)).transpose()
    if outname:
        np.savetxt(outname,result,fmt='%.2f')
    if plot:
        plt.plot(freqs,ratio,'.')
        plt.show()
    return result

if __name__ == '__main__':
    dists = ['20','30','40','50','60','70','80']
    for dist in dists:
        zfile = 'g%s.z' %dist
        rfile = 'g%s.r' %dist
        ex = cal_zhratio(zfile,rfile,np.arange(0.01,0.21,0.01))
        ex = ex.transpose()
        forward = np.loadtxt('zhr_forward')

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(forward[0],forward[1])
        ax.plot(ex[0],ex[1],'.',color='red')
        ax.set_xlabel('frequency(Hz)')
        ax.set_ylabel('Z/Hratio')
        ax.set_xlim(0.005,0.21)
        ax.set_ylim(1,1.6)
        ax.set_title('Z/H ratio plot of '+dist)
        plt.savefig(dist+'.png')
