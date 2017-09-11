from subprocess import Popen,PIPE
from glob import glob

def trans(sacfile):
    """
    sac file from seed
    rmean;rtrend;taper
    trans
    write staname.BH?
    """
    p = Popen(['sac'], stdin=PIPE, stdout=PIPE)
    temp = sacfile.split('.')
    net,sta,code,ch = temp[6],temp[7],temp[8],temp[9]
    respfile = '.'.join(['RESP',net,sta,code,ch])
    outfile = '.'.join(['SAC',net,sta,code,ch])
    f = [0.5,1,28,30]
    s = ""
    s += "r %s\n" % sacfile
    s += "rmean;rtrend;taper\n"
    s += "transfer from evalresp fname %s to none freq %f %f %f %f\n" % (respfile,f[0],f[1],    f[2],f[3])
    s += "w %s\n" % outfile
    s += "q\n"
    p.communicate(s.encode())

def get_three(file_list):
    result = []
    common = lambda s : '.'.join(s.split('.')[:-1])
    comm_parts = set(map(common,file_list))
    for cp in comm_parts:
        result.append(glob(cp+'*'))
    return result

def rotate(file_list2):
    for three in file_list2:
        two = filter(lambda s:True if s.split('.')[-1]!='BHZ' else False,three)
        two = list(two)
        p = Popen(['sac'], stdin=PIPE, stdout=PIPE)
        s = ""
        s += "r %s\n" % ' '.join(two)
        s += "rotate to gcp\n"
        comm = '.'.join(two[0].split('.')[:-1])
        s += "w %s.R %s.T\n" %(comm,comm)
        s += "q\n"
        p.communicate(s.encode())

def do_list(file_list):
    for sacfile in file_list:
        trans(sacfile)
