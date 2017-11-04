import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def do_mft(filename,wave,dist):
    chose_alpha = {1000:25,2000:50,4000:100,18000:200}
    command = ['sacmft96','-f']
    command.append(filename)
    for distmax in chose_alpha:
        if dist <= distmax:
            alpha = chose_alpha[distmax]
            break
    command.extend(['-a0',str(alpha)])
    command.extend(['-PMIN','5','-PMAX','120'])
    if wave == 'R':
        command.append('-R')
    elif wave == 'L':
        command.append('-L')
    subprocess.run(command,check=True)
    subprocess.run("awk '{print $5,$6,$10}' mft96.disp > temp",shell=True)
    df = pd.read_table('temp',sep='\s+',names=['per','vel','amp'])
    result = df.groupby('per').apply(lambda x:x['vel'][x['amp'].argmax()])
    #plt.plot(result)
    #plt.show()
    subprocess.run("rm mft* MFT* temp",shell=True)
    
    f = result.index
    f = 1.0/f
    c = result.values
    return np.vstack([f[::-1],c[::-1]])

if __name__ == '__main__':
    result = do_mft('g80.z','R',8895.59)
