from glob import glob
import numpy as np
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


def reject_outliers(prev_dir, out_dir):
    """
    statistically reject outliers
    :param prev_dir: directory of dispersion curve
    :param out_dir: directory to save selected dispersion,
        result[0] mean vel, result[1] std
    :return: None
    """
    thr_num = 10
    thr_len = 40
    PRANGE = (10, 80)
    verbose = True
    n = 1.5

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    prev_dirs = glob(prev_dir + '*')
    for dire in prev_dirs:
        pair = dire.split('/')[-1]
        results = []
        disp_files = glob(dire + '/*')
        if len(disp_files) < thr_num:
            continue
        for disp_file in disp_files:
            results.append(np.loadtxt(disp_file))
        results = np.array(results)
        mean = np.nanmean(results, axis=0)
        std = np.nanstd(results, axis=0)
        # quantity and quality(std) of data
        mask = (std / mean) < 0.03
        mask = mask & (np.sum(~np.isnan(results), axis=0) > thr_num)
        out_range = max_conti_true(mask, thr_len)
        if out_range[1] - out_range[0] == 0:
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
        sstd = np.nanstd(selected, axis=0)
        window(smean, out_range)
        window(sstd, out_range)
        # write
        np.savetxt(out_dir+pair, np.vstack([smean, sstd]))
        # verbose mode
        if verbose:
            fig, axes = plt.subplots(1, 2, sharex=True, sharey=True)
            fig.set_size_inches(10, 5)
            for i in range(len(results)):
                axes[0].plot(np.arange(*PRANGE), results[i], color='orange')
                axes[0].plot(np.arange(*PRANGE), selected[i], color='black')
            axes[1].plot(np.arange(*PRANGE), mean, color='red', label='raw')
            axes[1].plot(np.arange(*PRANGE), smean, color='blue', label='selected')
            axes[1].grid(color='grey', linestyle='dashed')
            axes[0].grid(color='grey', linestyle='dashed')
            axes[1].legend()
            axes[0].set_ylim(2.9, 4.3)
            axes[0].set_xlim(10, 80)
            plt.show(block=False)
            _ = input('>')
            plt.close()


if __name__ == '__main__':
    root = '/home/haosj/data/neTibet/'
    reject_outliers(root+'result/', root+'result2/')
