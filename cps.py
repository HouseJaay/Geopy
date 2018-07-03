import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from glob import glob


def do_mft(filename, wave, dist):
    """
    call sacmft96 from Computer Programs in Seismology
        calculate group velocity curve
    :param filename: sac file name
    :param wave: wave type ,'L' or 'R'
    :param dist: epicenter distance(km)
    :return: result[0] frequency, result[1] group velocity
    """
    chose_alpha = {1000: 25, 2000: 50, 4000: 100, 18000: 200}
    command = ['sacmft96', '-f', filename]
    for distmax in chose_alpha:
        if dist <= distmax:
            alpha = chose_alpha[distmax]
            break
    command.extend(['-a0', str(alpha)])
    command.extend(['-PMIN', '5', '-PMAX', '120'])
    if wave == 'R':
        command.append('-R')
    elif wave == 'L':
        command.append('-L')
    with open(os.devnull, 'w') as devnull:
        subprocess.run(command, check=True, stdout=devnull)
    subprocess.run("awk '{print $5,$6,$10}' mft96.disp > temp", shell=True)
    df = pd.read_table('temp', sep='\s+', names=['per', 'vel', 'amp'])
    result = df.groupby('per').apply(lambda x: x['vel'][x['amp'].argmax()])
    #plt.plot(result)
    #plt.show()
    subprocess.run("rm mft* MFT* temp", shell=True)
    f = result.index
    f = 1.0/f
    c = result.values
    return np.vstack([f[::-1], c[::-1]])


def litho_to_mod96(lat, lon, outname):
    """
    access LITHO1.0 and convert to mod96 format
    :param lat: latitude
    :param lon: longitude
    :param outname: mod96 model file name
    :return:
    """
    subprocess.run(
        "access_litho -p %d %d | awk '{print $1,$2,$3,$4,$5,$6}' > litho.temp"
        % (lat, lon), shell=True)
    model = np.loadtxt("litho.temp")
    model = model[::-1, :]
    model[:, :4] = model[:, :4]/1000
    model[:, 0] = model[:, 0] - model[0][0]
    # convert Qmu,Qkappa to Qp,Qs
    L = (4/3) * (model[:, 3]/model[:, 2])**2
    Qp = 1/L * model[:, 5] + 1/(1-L) * model[:, 4]
    Qs = model[:, 5]
    model[:, 4], model[:, 5] = Qp, Qs
    # change column sequence
    model[:, 1], model[:, 2], model[:, 3] = \
        model[:, 2].copy(), model[:, 3].copy(), model[:, 1].copy()
    # convert depth to layer thieckness
    cps_model = []
    if len(model) % 2 == 0:
        raise ValueError('unexpected model for lat:%d lon:%d' % (lat, lon))
    for i in range(0, len(model)-1, 2):
        layer_thickness = model[i+1, 0] - model[i, 0]
        if all(model[i, 1:] == model[i+1, 1:]):
            model[i, 0] = layer_thickness
            cps_model.append(model[i, :])
        else:
            raise ValueError('unexpected model for lat:%d lon:%d' % (lat, lon))
    model[-1, 0] = 100
    cps_model.append(model[-1, :])
    cps_model = np.array(cps_model)
    # add ETAP ETAS FREFP FREFS
    comp = np.zeros([len(cps_model), 4])
    comp[:, 2:] = 1
    cps_model = np.concatenate((cps_model, comp), axis=1)
    
    header = "%s\nlitho1.0\n0" % outname
    np.savetxt("model.temp", cps_model, header=header, comments="", fmt='%8.4f')
    subprocess.run("cat model.temp | mkmod96", shell=True)
    subprocess.run("rm litho.temp model.temp", shell=True)
    return cps_model


def forward_rayleigh(modelname):
    """
    using given model, forward compute surface wave phase velocity
        group velocity and ZHratio
    :param modelname: mod96 format model
    :return: result[0] frequency, result[1] phase velocity,
        result[2] group velocity, result[3] ZHratio
    """
    subprocess.run(
        "sprep96 -M %s -HS 5 -HR 0 -DT 0.5 -NPTS 2048 -R -L -NMOD 1" % modelname,
        shell=True)
    subprocess.run('sdisp96', shell=True)
    subprocess.run('sregn96', shell=True)
    subprocess.run('slegn96', shell=True)
    subprocess.run(
        "sdpegn96 -R -C -U -PER -YMIN 2 -YMAX 5 -XMIN 1 -XMAX 80 -ASC",
        shell=True)
    subprocess.run("awk '{print $4,$5,$6,$9}' SREGN.ASC > temp", shell=True)
    result = np.loadtxt('temp', skiprows=1)
    subprocess.run(
        'rm sdisp* temp slegn96.egn sregn96.egn SREGN*', shell=True)
    
    result = result.transpose()
    result[3] = 1 / result[3]
    return result[:, ::-1]  # frequency ascending order


def forward_love(modelname):
    """
    :param modelname:
    :return: result[0] frequency, result[1] phase velocity,
        result[2] group velocity
    """
    subprocess.run(
        "sprep96 -M %s -HS 5 -HR 0 -DT 0.5 -NPTS 2048 -R -L -NMOD 1" % modelname,
        shell=True)
    subprocess.run('sdisp96', shell=True)
    subprocess.run('sregn96', shell=True)
    subprocess.run('slegn96', shell=True)
    subprocess.run(
        "sdpegn96 -L -C -U -PER -YMIN 2 -YMAX 5 -XMIN 1 -XMAX 80 -ASC",
        shell=True)
    subprocess.run("awk '{print $4,$5,$6}' SLEGN.ASC > temp", shell=True)
    result = np.loadtxt('temp', skiprows=1)
    subprocess.run(
        'rm sdisp* temp slegn96.egn sregn96.egn SLEGN*', shell=True)

    result = result.transpose()
    return result[:, ::-1]  # frequency ascending order


def compute2d(lat, lon, mark, outname):
    """
    compute 2d data
    :param lat: latitude range
    :param lon: longitude range
    :param mark: 1 phase vel, 2 group vel, 3 ZHratio
    :param outname: write filename
    :return:
    """
    freqs = [0.1, 0.05, 0.04, 0.025, 0.02, 0.0125]
    out = ""
    for la in range(*lat):
        for lo in range(*lon):
            litho_to_mod96(la, lo, 'temp')
            forward = forward_rayleigh('temp')
            results = np.interp(freqs, forward[0], forward[mark])
            out += str(lo) + ' ' + str(la) + ' '
            out += ' '.join(map(lambda x: '%.4s' % x, results)) + '\n'
    with open(outname, 'w') as f:
        f.write(out)


def plot_mod96(filename):
    fig, ax = plt.subplots()
    with open(filename, 'r') as f:
        models = f.readlines()
    vp, vs, h = [], [], []
    for layer in models[12:]:
        vp.append(float(layer.split()[1]))
        h.append(float(layer.split()[0]))
        vs.append(float(layer.split()[2]))
    if h[-1] == 0:
        h[-1] = sum(h)*0.1  # infinite half space
    vpd, vsd, depth = [], [], []
    for i in range(len(h)):
        if not depth:
            depth.append(0)
        else:
            depth.append(depth[-1])
        vpd.append(vp[i])
        vsd.append(vs[i])
        depth.append(depth[-1]+h[i])
        vpd.append(vp[i])
        vsd.append(vs[i])
    plt.title(models[1][:-1], y=1.08)
    ax.set_xlim([min(vs)-1, max(vp)+1])
    ax.set_xlabel("velocity (km/s)")
    ax.set_ylim(0, sum(h))
    ax.set_ylabel("depth (km)")
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax.plot(vsd, depth, 'red', label='vs')
    ax.plot(vpd, depth, 'blue', label='vp')
    plt.legend()
    plt.show()


class Inv:
    @staticmethod
    def write_surf96(filename, wave, type, flag, mode, peroid_arr, value_arr, err_arr):
        """
        write to surf96 format dispersion file
        :param :refer to cps doc
        :param filename: filename to write
        """
        f = open(filename, 'w')
        for i in range(len(peroid_arr)):
            f.write("SURF96 %s %s %s %d %f %f %f\n" % (
                wave, type, flag, mode, peroid_arr[i], value_arr[i], err_arr[i]
            ))
        f.close()

    @staticmethod
    def read_disp(filename, peroid):
        disp = np.loadtxt(filename)
        p_disp = np.c_[peroid.T, disp.T][~np.isnan(disp[0])]
        Inv.write_surf96('disp.d', 'R', 'C', 'X', 0, p_disp[:, 0],
                         p_disp[:, 1], p_disp[:, 2])

    @staticmethod
    def do_inv_netibet(datadir, outdir):
        iter_num = 5

        disps = glob(datadir+'*')
        for disp in disps:
            name = disp.split('/')[-1]
            Inv.read_disp(disp, np.arange(10, 80))
            # TODO write model
            # TODO change parameter
            subprocess.run("surf96"+" 1 2 6"*iter_num, shell=True)
            subprocess.run("surf96 28 %s" % outdir+name, shell=True)
            subprocess.run("surf96 39", shell=True)
            # TODO monitor output

    @staticmethod
    def conv_to_plain_vsmodel(datadir, outdir):
        models = glob(datadir+'*')
        for model in models:
            name = model.split('/')[-1]
            with open(model, 'r') as f:
                lines = f.readlines()[12:]
            vs = ""
            err = ""
            for line in lines:
                vs += line.split()[2] + ' '
                err += '0.05 '  # TODO change it
            with open(outdir+name, 'w') as f:
                f.write(vs + '\n')
                f.write(err + '\n')


if __name__ == '__main__':
    model = litho_to_mod96(30, 108, "./testdata/ZHratio/test.d")
    groupv = do_mft('./testdata/ZHratio/g30.r', 'R', 3335.8)
    forward = forward_rayleigh("./testdata/ZHratio/test.d")
    compute2d((32, 42), (96, 108), 1, 'litho1.0_phasevel')
