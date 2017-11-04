from scipy import signal
import obspy
import matplotlib.pyplot as plt 
import numpy as np
from Geopy import metadata as mt
from Geopy import cps
from obspy.signal.filter import envelope

def plot_progress(Z,Z_env,H,H_env,freqs,Z_win,H_win):
    fig,axes = plt.subplots(len(freqs)//2+1,2)
    for i in range(len(freqs)):
        if i%2 == 0:
            j = i//2
            axes[j,0].plot(Z[i])
            axes[j,0].plot(Z_env[i])
            axes[j,1].plot(H[i])
            axes[j,1].plot(H_env[i])
            axes[j,0].set_ylabel('freq'+str(freqs[i]))
            axes[j,0].plot(Z_win[i]*Z[i].max())
            axes[j,1].plot(H_win[i]*H[i].max())
    plt.show()

def group_vel_win(filename,stats,freqs,n):
    dist = stats.sac['dist']
    delta = stats.delta
    timediff = stats.sac['b'] - stats.sac['o']
    npts = stats.npts
    gvdisp = cps.do_mft(filename,'R',dist)
    gvinterp = np.interp(freqs,gvdisp[0],gvdisp[1])
    pers = 1/freqs
    left = (dist/(gvinterp+1) - timediff)/delta
    left = left.astype('int')
    right = (dist/(gvinterp-1) -timediff)/delta
    right = right.astype('int')
    def cut(bound):
        win = np.zeros(npts)
        l = 0 if bound[0]<0 else bound[0]
        r = npts-1 if bound[1]>npts-1 else bound[1]
        win[l:r] = 1
        return win
    return list(map(cut,zip(left,right)))
    
def cal_zhratio(zfile,rfile,freqs,bpwidth=0.002,outname=None,plot=None):
    try:
        st = obspy.read(zfile)
        st += obspy.read(rfile)
    except Exception:
        print("file error %s %s" %(zfile,rfile))
        return 0
    delta = st[0].stats.delta
    dist = st[0].stats.sac['dist']
    def_b = lambda f:signal.firwin(1001,[f-bpwidth,f+bpwidth],
    window=('kaiser',9),nyq=1/delta/2,pass_zero=False)
    B = list(map(def_b,freqs))
    Z_win = group_vel_win(zfile,st[0].stats,freqs,3)
    H_win = group_vel_win(rfile,st[0].stats,freqs,3)
    Z = list(map(lambda b:signal.lfilter(b,1,st[0].data),B))
    H = list(map(lambda b:signal.lfilter(b,1,st[1].data),B))
    Z_env = list(map(envelope,Z))
    Z_env = list(map(lambda x:x[0]*x[1],zip(Z_env,Z_win)))
    H_env = list(map(envelope,H))
    H_env = list(map(lambda x:x[0]*x[1],zip(H_env,H_win)))

    ratio = np.array(list(map(lambda e:max(e[0])/max(e[1]),
    zip(Z_env,H_env))))

    result = np.vstack((freqs,ratio)).transpose()
    if outname:
        np.savetxt(outname,result,fmt='%.2f')
    if plot==1:
        plt.plot(freqs,ratio,'.')
        plt.show()
    if plot==2:
        plot_progress(Z,Z_env,H,H_env,freqs,Z_win,H_win)
    return result

if __name__ == '__main__':
    if False:
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
    if True:
        freqs = np.arange(0.01,0.21,0.01)
        cal_zhratio('g20.z','g20.r',freqs,plot=2)
