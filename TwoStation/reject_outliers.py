from glob import glob
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import os


def max_conti_true(arr, thr):
    """
    find longest continuous true subarray which length > thr
    :param arr: boolean array
    :param thr: minimum length, suppose thr > len(arr)/2
    :return: subarray index boundary [l, r]
        if no such subarray, return [0, 0]
    """
    j = 0
    while j < len(arr):
        while j < len(arr) and (not arr[j]):
            j += 1
        k = j
        while j < len(arr) and arr[j]:
            j += 1
        if j - k > thr:
            return [k, j]
    return [0, 0]


def window(arr, win):
    arr[:win[0]] = np.NaN
    arr[win[1]:] = np.NaN


def reject_outliers(prev_dir, out_dir, thr_num=5, thr_len=40, thr_std=0.03,
                    x_value=np.arange(10, 80), n=1.5, verbose=False):
    """
    statistically reject outliers
    :param prev_dir: directory of dispersion curve
    :param out_dir: directory to save selected dispersion,
        result[0] mean vel, result[1] std
    :param thr_num: threshold quantities of result
    :param thr_len: threshold quantities of continuous data point
    :param thr_std: threshold of std/mean
    :param x_value: frequency of peroids corresponding to result data, used to plot
    :param n: control outlier threshold
    :param verbose: display detailed plot
    :return: None
    """
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    prev_dirs = glob(prev_dir + '*')
    for dire in prev_dirs:
        print(dire)
        pair = dire.split('/')[-1]
        results = []
        disp_files = glob(dire + '/*')
        if len(disp_files) < thr_num:
            print('1')
            continue
        for disp_file in disp_files:
            results.append(np.loadtxt(disp_file))
        results = np.array(results)
        mean = np.nanmean(results, axis=0)
        std = np.nanstd(results, axis=0)
        # quantity and quality(std) of data
        print(std/mean)
        mask = (std / mean) < thr_std
        mask = mask & (np.sum(~np.isnan(results), axis=0) > thr_len)
        out_range = max_conti_true(mask, thr_len)
        if out_range[1] - out_range[0] == 0:
            print('2')
            continue
        # reject outliers
        selected = np.zeros(results.shape)
        selected[:, :] = np.NaN
        for i in range(len(disp_files)):
            mask = abs(results[i] - mean) < n * std
            l, r = max_conti_true(mask, thr_len)
            selected[i, l:r] = results[i, l:r]
        print(pair)
        # calculate mean and std of selected data, and window it
        smean = np.nanmean(selected, axis=0)
        # sstd = np.nanstd(selected, axis=0)
        sstd = stats.sem(selected, axis=0, nan_policy='omit')
        window(smean, out_range)
        window(sstd, out_range)
        # write
        np.savetxt(out_dir+pair, np.vstack([smean, sstd]))
        # verbose mode
        if verbose:
            ylim = (1.0, 1.7)
            # ylim = (2.9, 4.3)
            fig, axes = plt.subplots(1, 2, sharex=True, sharey=True)
            fig.set_size_inches(10, 5)
            for i in range(len(results)):
                axes[0].plot(x_value, results[i], color='orange')
                axes[0].plot(x_value, selected[i], color='black')
            axes[1].plot(x_value, mean, color='red', label='raw')
            axes[1].plot(x_value, smean, color='blue', label='selected')
            axes[1].grid(color='grey', linestyle='dashed')
            axes[0].grid(color='grey', linestyle='dashed')
            axes[1].legend()
            axes[0].set_ylim(*ylim)
            axes[0].set_xlim(min(x_value), max(x_value))
            plt.show(block=False)
            _ = input('>')
            plt.close()


if __name__ == '__main__':
    root = '/home/haosj/data/neTibet/'
    reject_outliers(root+'result/', root+'result_new/')
