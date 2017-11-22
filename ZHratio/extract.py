from scipy import signal
import obspy
import matplotlib.pyplot as plt 
import numpy as np
from Geopy import metadata as mt
from Geopy import cps
from obspy.signal.filter import envelope

def plot_progress(zdata,hdata,Z,Z_env,H,H_env,freqs,Z_win,H_win):
    fig,axes = plt.subplots(len(freqs)//2+2,2)
    axes[0,0].plot(zdata)
    axes[0,1].plot(hdata)
    for i in range(len(freqs)):
        if i%2 == 0:
            j = i//2+1
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
    left = (dist/(gvinterp+n) - timediff)/delta
    left = left.astype('int')
    right = (dist/(gvinterp-n) -timediff)/delta
    right = right.astype('int')
    def cut(bound):
        win = np.zeros(npts)
        l = 0 if bound[0]<0 else bound[0]
        r = npts-1 if bound[1]>npts-1 else bound[1]
        win[l:r] = 1
        return win
    return list(map(cut,zip(left,right)))

def sel_evt(st):
    def helper(tr):
        if tr.stats.sac['mag'] > 6.0 and tr.stats.sac['evdp'] <= 20 and\
        tr.stats.sac['dist'] > 1500 and tr.stats.sac['dist'] < 8000:
            start = tr.stats.starttime
            return "X2.%s.%d.%03d.%02d.%02d.%02d.00" %(
            tr.stats.station,start.year,start.julday,start.hour,
            start.minute,start.second)
        else:
            return None
    return list(filter(lambda x:True if x else False,map(helper,st)))
 
def check_corr(Z,H,freqs,delta):
    """
    vertical component phase advance 90 dgrees
    """
    norm = list(map(lambda x:(np.dot(x[0],x[0])*np.dot(x[1],x[1]))**0.5
    ,zip(Z,H)))
    COR = list(map(lambda x:signal.correlate(x[0],x[1])/x[2],zip(Z,H,norm)))
    fig,axes = plt.subplots(len(freqs),1,sharex=True)
    for i,cor in enumerate(COR):
        start = int((len(cor)-1)/2)
        t = 1/freqs[i]
        n = int(t/4/delta)
        end = start + 4*int(t/delta)
        cor = cor[start:end]
        marker = np.zeros(len(cor))
        marker[n] = 0.8
        print(cor[n])
        axes[i].plot(cor,color='b')
        axes[i].plot(marker,color='r')
        axes[i].set_ylabel('%f' %freqs[i])
    plt.show()

def check_corr2(Z,H,freqs,delta):
    norm = list(map(lambda x:(np.dot(x[0],x[0])*np.dot(x[1],x[1]))**0.5
    ,zip(Z,H)))
    COR = list(map(lambda x:signal.correlate(x[0],x[1])/x[2],zip(Z,H,norm)))
    def helper(x):
        cor,freq = x[0],x[1]
        n = int((len(cor)-1)/2) + int(1/freq/4/delta)
        return cor[n]
    return np.array(list(map(helper,zip(COR,freqs))))

def cal_zhratio(zfile,rfile,freqs,bpwidth=0.002,outname=None,plot=None,
threshold=0.8):
    try:
        st = obspy.read(zfile)
        st += obspy.read(rfile)
    except Exception:
        print("file error %s %s" %(zfile,rfile))
        return 0
    st[0].data = signal.detrend(st[0].data)
    st[1].data = signal.detrend(st[1].data)
    delta = st[0].stats.delta
    dist = st[0].stats.sac['dist']
    def_b = lambda f:signal.firwin(1001,[f-bpwidth,f+bpwidth],
    window=('kaiser',9),nyq=1/delta/2,pass_zero=False)
    B = list(map(def_b,freqs))
    Z_win = group_vel_win(zfile,st[0].stats,freqs,1)
    H_win = group_vel_win(rfile,st[0].stats,freqs,1)
    Z = list(map(lambda b:signal.lfilter(b,1,st[0].data),B))
    Z = list(map(lambda x:x[0]*x[1],zip(Z_win,Z)))
    H = list(map(lambda b:signal.lfilter(b,1,st[1].data),B))
    H = list(map(lambda x:x[0]*x[1],zip(H_win,H)))
    Z_env = list(map(envelope,Z))
    H_env = list(map(envelope,H))
    
    cor_eff = check_corr2(Z,H,freqs,delta)
    if len(cor_eff[cor_eff>threshold]) < 3:
        return False

    ratio = np.array(list(map(lambda e:max(e[0])/max(e[1]),
    zip(Z_env,H_env))))
    result = np.vstack((freqs,ratio)).transpose()
    if outname:
        np.savetxt(outname,result,fmt='%.2f')
    if plot==1 or plot==2:
        plt.plot(freqs,ratio,'.')
        plt.show()
    if plot==2:
        plot_progress(st[0].data,st[1].data,
        Z,Z_env,H,H_env,freqs,Z_win,H_win)
    return np.array(list(map(lambda x:x[0] if x[1]>threshold else None,
    zip(ratio,cor_eff))),dtype=np.float)

if __name__ == '__main__':
    if False: # plot synthetic test
        dists = ['20','30','40','50','60','70','80']
        rows,cols = len(dists)//2+1,2
        fig = plt.figure()
        for i,dist in enumerate(dists):
            ax = fig.add_subplot(rows,cols,i+1)
            zfile = 'g%s.z' %dist
            rfile = 'g%s.r' %dist
            ex = cal_zhratio(zfile,rfile,np.arange(0.01,0.21,0.01))
            ex = ex.transpose()
            forward = np.loadtxt('zhr_forward')

            ax.plot(forward[0],forward[1])
            ax.plot(ex[0],ex[1],'.',color='red')
            ax.set_ylabel('Z/Hratio '+dist)
            ax.set_xlim(0.005,0.21)
            ax.set_ylim(1,1.6)
            ax.set_xlabel('frequency(Hz)')
        plt.show()
        #plt.savefig('synthetic.png')
    if False:
        freqs = np.arange(0.01,0.21,0.01)
        cal_zhratio('g80.z','g80.r',freqs,plot=2)
    if True:

        st = obspy.read('/home/haosj/data/tibet/ped/15639/*.Z')
        commons = sel_evt(st)
        freqs = np.arange(0.01,0.11,0.01)
        results = []
        for common in commons:
            temp = cal_zhratio(common+'.Z',common+'.R',freqs,plot=0,
            threshold=0.8)
            if isinstance(temp,np.ndarray):
                results.append(temp)
        results = np.array(results)
        for i in range(len(freqs)):
            if np.count_nonzero(~np.isnan(results[:,i]))<3:
                results[:,i] = np.nan
        mean = np.nanmean(results,axis=0)
        std = np.nanstd(results,axis=0)
        fig,ax = plt.subplots()
        ax.errorbar(freqs,mean,yerr=std)
        plt.plot()
        plt.show()
