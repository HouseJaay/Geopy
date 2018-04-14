from subprocess import Popen,PIPE
from glob import glob
import obspy
import os
from os.path import isfile

os.putenv("SAC_DISPLAY_COPYRIGHT",'0')

def trans(filename,respname,writename,f,d):
    p = Popen(['sac'], stdin=PIPE, stdout=PIPE)
    s = ""
    s += "r %s\n" % filename
    s += "decimate %d;decimate %d\n" % (d[0],d[1])
    s += "rmean;rtrend\n"
    s += "transfer from polezero subtype %s to none freq %f %f %f %f\n" % (
    respname,f[0],f[1],f[2],f[3])
    s += "w %s\n" % writename
    s += "q\n"
    r = p.communicate(s.encode())
    #print(r)

def do_trans(sacdir,respdir,writedir):
    saclist = glob(sacdir+'*.SAC')
    print('processing %s' %sacdir)
    for sacpath in saclist:
        sacfile = sacpath.split('/')[-1]
        sta,ch = sacfile.split('.')[1],sacfile.split('.')[-2]
        respname = respdir + "SAC_PZs_X2_%s_%s_00" %(sta,ch)
        writename = writedir + '.'.join(sacfile.split('.')[:-2]) +\
        '.' + ch[-1]
        trans(sacpath,respname,writename,(0.008,0.012,3,4),(5,2))

def do_rotate(peddir):
    zfiles = glob(peddir+'*.Z')
    for zfile in zfiles:
        common = '.'.join(zfile.split('.')[:-1])
        efile,nfile = common + '.E',common + '.N'
        rfile,tfile = common + '.R',common + '.T'
        if isfile(efile) and isfile(nfile):
            p = Popen(['sac'],stdin=PIPE,stdout=PIPE)
            s = ""
            s += "r %s\n" %efile
            s += "ch cmpinc 90 cmpaz 90;wh\n"
            s += "r %s\n" %nfile
            s += "ch cmpinc 90 cmpaz 0;wh\n"
            s += "r %s %s\n" %(efile,nfile)
            s += "rotate to gcp\n"
            s += "w %s %s\n" %(rfile,tfile)
            s += "q\n"
            r = p.communicate(s.encode())
            #print(r)
        else:
            print('file error %s' %zfile)
    
if __name__ == '__main__':
    root = '/home/haosj/data/tibet/'
    dirs = os.listdir(root+'ordos/')
    for sta in dirs:
        writedir = root + 'ped/' + sta + '/'
        os.mkdir(writedir)
        do_trans(root+'ordos/'+sta+'/',root+'RESP/RESP/',writedir)
        do_rotate(writedir)
