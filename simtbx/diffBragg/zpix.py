import numpy as np
import glob
from iotbx.reflection_file_reader import any_reflection_file
import sys
import os

ftruth = any_reflection_file(sys.argv[1]).as_miller_arrays()[0]
if not ftruth.is_xray_amplitude_array():
    ftruth = ftruth.as_amplitude_array()
res = ftruth.d_spacings()
res_map = {h: d for h, d in zip(res.indices(), res.data())}
datadir = sys.argv[2]
asu_map = np.load(os.path.join(datadir, "f_asu_map.npy"), allow_pickle=True)[()]

measured_res = np.sort([res_map[h] for h in asu_map.values()])
nbins = 9
res_bins = [s[0] for s in np.array_split(measured_res, nbins)] + [measured_res[-1]]
res_bin_cent = [(a+b)/2. for a,b in zip(res_bins, res_bins[1:])]

Zscore_dirs = glob.glob("%s/Z/rank*_Zscore" % datadir)
import time

all_res =[]
all_sigZ = []
all_iter = []

maxiters = 500
img_data = np.zeros((nbins, maxiters))
img_data_count = np.zeros((nbins, maxiters))
for i_dir, dirname in enumerate(Zscore_dirs):
    t = time.time()
    fnames = glob.glob("%s/*npy" % dirname)
    print("Found %d fnbames in dir %s" % (len(fnames), dirname))
    iters = [int(os.path.basename(f).split("_")[1].split("iter")[1]) for f in fnames]
    iters, fnames = zip(*sorted(zip(iters, fnames)))

    for i_f, f in enumerate(fnames):
        rank_data = np.load(f, allow_pickle=True)[()]

        for shot_data in rank_data:
            i_fcells, sigZ = map(np.array, zip(*shot_data))
            res = np.array([res_map[asu_map[i_fcell]] for i_fcell in i_fcells])
            #all_res += res
            #all_sigZ += sigZ
            #all_iter += [iters[i_f]] * len(res)
            digs = np.digitize(res, res_bins)
            for i_bin in range(nbins):
                sel = digs == (i_bin+1)
                img_data[i_bin, iters[i_f]] += np.sum(sigZ[sel])
                img_data_count[i_bin, iters[i_f]] += np.sum(sel)

    t = time.time()-t
    print("took %f sec" % t)
    #if i_dir == 10:
    #    break

from pylab import *
D = np.nan_to_num(img_data / img_data_count)
all_iter = np.arange(D.shape[1])
for i_bin in range(len(D)):
    d = D[i_bin]
    sel = d > 0
    iters = all_iter[sel]
    a, b = res_bins[i_bin], res_bins[i_bin+1]
    if i_bin==0:
        label="%.1f - %.1f $\AA$" % (a,b)
    else:
        label = "%.1f - %.1f" % (a, b)
    plot(iters, d[sel], marker='.',label=label)
legend(ncol=3, markerscale=1.5)
ax = gca()
ax.tick_params(labelsize=11)
xlabel("stage 2 L-BFGS iteration", fontsize=14)
ylabel("average Z-score stdev.", fontsize=14)
grid(1, alpha=0.5)
title(sys.argv[2])
show()

