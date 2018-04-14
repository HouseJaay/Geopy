from scipy import signal
from scipy.fftpack import hilbert
import obspy
import matplotlib.pyplot as plt
import numpy as np
from Geopy import metadata as mt
from Geopy import cps
from obspy.signal.filter import envelope
import pandas
import os

def plot_progress(zdata, hdata, Z, Z_env, H, H_env, freqs, Win):
    fig, axes = plt.subplots(len(range(0, len(freqs), 2))+1, sharex=True)
    zcolor, hcolor = 'red', 'blue'
    axes[0].plot(zdata, color=zcolor)
    axes[0].plot(hdata, color=hcolor)
    for i in range(0, len(freqs), 2):
        j = i//2+1
        axes[j].plot(Z[i], color=zcolor)
        axes[j].plot(Z_env[i], color=zcolor)
        axes[j].plot(H[i], color=hcolor)
        axes[j].plot(H_env[i], color=hcolor)
        axes[j].set_ylabel('freq_'+str(freqs[i]))
        axes[j].plot(Win[i]*Z[i].max())
    plt.show()

def group_vel_win(filename, stats, freqs, n):
    dist = stats.sac['dist']
    delta = stats.delta
    timediff = stats.sac['b'] - stats.sac['o']
    npts = stats.npts
    gvdisp = cps.do_mft(filename, 'R', dist)
    gvinterp = np.interp(freqs, gvdisp[0], gvdisp[1])
    pers = 1/freqs
    left = (dist/(gvinterp+n) - timediff)/delta
    left = left.astype('int')
    right = (dist/(gvinterp-n) -timediff)/delta
    right = right.astype('int')
    def cut(bound):
        win = np.zeros(npts)
        l = 0 if bound[0] < 0 else bound[0]
        r = npts-1 if bound[1] > npts-1 else bound[1]
        win[l:r] = 1
        return win
    return list(map(cut, zip(left, right)))

def sel_evt(st):
    def helper(tr):
        if tr.stats.sac['mag'] > 5 and tr.stats.sac['evdp'] <= 40 and\
        tr.stats.sac['dist'] > 1500 and tr.stats.sac['dist'] < 8000:
            start = tr.stats.starttime
            return "X2.%s.%d.%03d.%02d.%02d.%02d.00" %(
            tr.stats.station,start.year,start.julday,start.hour,
            start.minute,start.second)
        else:
            return None
    return list(filter(lambda x:True if x else False,map(helper,st)))
 
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
    st[0].data = signal.detrend(st[0].data) #hilbert transform z component
    st[0].data = hilbert(st[0].data)
    st[1].data = signal.detrend(st[1].data)
    delta = st[0].stats.delta
    dist = st[0].stats.sac['dist']
    def_b = lambda f:signal.firwin(1001,[f-bpwidth,f+bpwidth],
    window=('kaiser',9),nyq=1/delta/2,pass_zero=False)
    B = list(map(def_b,freqs))
    Win = group_vel_win(zfile,st[0].stats,freqs,1)
    Z = list(map(lambda b:signal.lfilter(b,1,st[0].data),B))
    Z = list(map(lambda x:x[0]*x[1],zip(Win,Z)))
    H = list(map(lambda b:signal.lfilter(b,1,st[1].data),B))
    H = list(map(lambda x:x[0]*x[1],zip(Win,H)))
    Z_env = list(map(envelope,Z))
    H_env = list(map(envelope,H))
    
    cor_eff = np.array(list(map(lambda x:np.corrcoef(
    x[0],x[1])[1,0],zip(H,Z))))
    if len(cor_eff[cor_eff>threshold]) < 3:
        return False

    ratio = np.array(list(map(lambda e:max(e[0])/max(e[1]),
    zip(Z_env,H_env))))
    result = np.vstack((freqs,ratio)).transpose()
    if outname:
        np.savetxt(outname,result,fmt='%.2f')
    if plot==1:
        plt.plot(freqs,ratio,'.')
        plt.show()
    if plot==2:
        plot_progress(st[0].data,st[1].data,
        Z,Z_env,H,H_env,freqs,Win)
    return np.array(list(map(lambda x:x[0] if x[1]>threshold else None,
    zip(ratio,cor_eff))),dtype=np.float)

def do_single_sta(sta_dir,freqs):
    st = obspy.read(sta_dir + '*.Z',headonly=True)
    commons = sel_evt(st)
    print(len(commons))
    baz,results = [],[]
    for common in commons:
        fname = sta_dir + common
        try:
            temp = cal_zhratio(fname+'.Z',fname+'.R',freqs,plot=0,threshold=0.8)
        except Exception as e:
            print(fname+':')
            print(e)
            temp = False
        if isinstance(temp,np.ndarray):
            results.append(temp)
            st = obspy.read(fname+'.Z',headonly=True)
            baz.append(st[0].stats.sac['baz'])
    ratio_baz = sorted(zip(results,baz), key=lambda x:x[1], reverse=True)
    deg_bin = 30
    cent_bazs,mean_results,std_results = [],[],[]
    for freq_rbound in range(deg_bin,360,deg_bin):
        cur_results = []
        while ratio_baz and ratio_baz[-1][1] < freq_rbound:
            cur_results.append(ratio_baz.pop()[0])
        if len(cur_results) > 1:
            cent_bazs.append(freq_rbound - deg_bin//2)
            mean_results.append(np.nanmean(cur_results, axis=0))
            std_results.append(np.nanstd(cur_results, axis=0))
    return results,baz,cent_bazs,mean_results,std_results
if __name__ == '__main__':
    if False: # plot synthetic test
        dists = ['20','30','40','50','60','70','80']
        rows,cols = len(dists)//2+1,2
        fig = plt.figure()
        freqs = np.arange(0.01,0.21,0.01)
        for i,dist in enumerate(dists):
            ax = fig.add_subplot(rows,cols,i+1)
            zfile = 'g%s.z' %dist
            rfile = 'g%s.r' %dist
            ex = cal_zhratio(zfile,rfile,freqs,threshold=0.8)
            forward = np.loadtxt('zhr_forward')

            ax.plot(forward[0],forward[1])
            ax.plot(freqs,ex,'.',color='red')
            ax.set_ylabel('Z/Hratio '+dist)
            ax.set_xlim(0.005,0.21)
            ax.set_ylim(1,1.6)
            ax.set_xlabel('frequency(Hz)')
        plt.show()
        #plt.savefig('synthetic.png')
    if False:
        freqs = np.arange(0.01,0.21,0.01)
        cal_zhratio('g30_cut.z','g30_cut.r',freqs,plot=2)
    if False:
        root = '/home/haosj/data/tibet/'
        sta = pandas.read_table(root+'metadata/ordos_sta.lst',
        names=['name','lat','lon','height'],sep='\s+',
        dtype={'name':str,'lat':np.float,'lon':np.float64})
        Mean,Std = [],[]
        mask = [True for i in range(len(sta))]
        for i in range(sta.shape[0]):
            filename = '%sped/%s/*.Z' %(root,sta['name'][i])
            st = obspy.read(filename,headonly=True)
            commons = sel_evt(st)
            freqs = np.arange(0.01,0.11,0.01)
            results = []
            for common in commons:
                fname = '%sped/%s/%s' %(root,sta['name'][i],common)
                temp = cal_zhratio(fname+'.Z',fname+'.R',freqs,plot=0,
                threshold=0.8)
                if isinstance(temp,np.ndarray):
                    results.append(temp)
            if len(results)>1:
                results = np.array(results)
                for i in range(len(freqs)):
                    if np.count_nonzero(~np.isnan(results[:,i]))<3:
                        results[:,i] = np.nan
                Mean.append(np.nanmean(results,axis=0))
                Std.append(np.nanstd(results,axis=0))
            else:
                mask[i] = False
        sta = sta[mask]
        Mean = np.array(Mean)
        Std = np.array(Std)
        lat = np.array(sta['lat'])
        lon = np.array(sta['lon'])
        zh = Mean[:,5]
    if True:
        freqs = np.arange(0.02,0.11,0.01)
        out = '/home/haosj/seis/zhratio/anis/'
        #stas = os.listdir('/home/haosj/data/tibet/ped/')
        stas = ['64046','64050','64053']
        for sta in stas:
            results,baz,cent_bazs,results_mean,results_std = do_single_sta('/home/haosj/data/tibet/ped/%s/' %sta,freqs)
            mean,std = np.array(results_mean),np.array(results_std)
            results = np.array(results)
            outdir = out+sta+'/'
            os.makedirs(outdir,exist_ok=True)
            for i,freq in enumerate(freqs):    
                fig, ax = plt.subplots(1)
                ax.plot(cent_bazs, mean[:,i], '.')
                ax.errorbar(cent_bazs, mean[:,i], yerr=std[:,i], fmt="none")
                ax.set_xlim(0,360)
                fmean = 'mean_' + str(freq) + '.png'
                plt.savefig(outdir+fmean)

                fig, ax = plt.subplots(1)
                ax.plot(baz,results[:,i],'.')
                ax.set_xlim(0,360)
                fscat = 'all_' + str(freq) +'.png'
                plt.savefig(outdir+fscat)


