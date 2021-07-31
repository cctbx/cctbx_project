import os
from cxid9114 import utils
from iotbx.reflection_file_reader import any_reflection_file
import numpy as np
import sys
from dials.array_family import flex
from cctbx import miller

ftruth = any_reflection_file(sys.argv[1]).as_miller_arrays()[0]
if not ftruth.is_xray_amplitude_array():
    ftruth = ftruth.as_amplitude_array()

datadir = sys.argv[2]  # = np.load(sys.argv[1])

asu_map = np.load(os.path.join(datadir, "f_asu_map.npy"), allow_pickle=True)[()]

fcell_pos, asu_indices = zip(*asu_map.items())
Fidx = flex.miller_index(asu_indices)

import glob

data_files = glob.glob(os.path.join(datadir, "_fcell_trial0_iter*.npz"))
n_files = len(data_files)
print("Found %d datafiles" % n_files)

data_files = sorted(data_files, key=lambda x: \
    int(os.path.basename(x).split("_iter")[1].split(".npz")[0]))

STRIDE = 40
indices = [0] + list(range(1, n_files - 2, STRIDE)) + [n_files - 1]

ftruth_map = None
mtz_names = []
all_data = []
for i_iter in indices:

    Fdata = flex.double(np.load(data_files[i_iter], allow_pickle=True)['fvals'])

    mset = ftruth.miller_set(flex.miller_index(Fidx), True)

    fobs = miller.array(mset, flex.double(Fdata))

    fobs_map = {h: d for h, d in zip(fobs.indices(), fobs.data())}

    if ftruth_map is None:
        hcommon = set(Fidx).intersection(ftruth.indices())
        ftruth_map = {h: d for h, d in zip(ftruth.indices(), ftruth.data())}
        hcommon_flex = flex.miller_index(list(hcommon))
        ftruth = miller.array(ftruth.miller_set(hcommon_flex, True), flex.double([ftruth_map[h] for h in hcommon_flex]))

    fobs = miller.array(fobs.miller_set(hcommon_flex, True),
                        flex.double([fobs_map[h] for h in hcommon_flex]))

    all_data.append(fobs.data().as_numpy_array())
    mtz_name = data_files[i_iter]+ ".mtz"
    mtz_names.append(mtz_names)
    fobs.as_mtz_dataset(column_root_label="F").mtz_object().write(mtz_name)
    out = utils.compute_r_factor(ftruth, fobs, verbose=False, is_flex=True, sort_flex=True)
    print(out)


from pylab import *
a = all_data[0]
b = all_data[-1]
print("Number of 0-valued Fhkl at start %f" % sum(a==0))
print("Number of 0-valued Fhkl at finish %f" % sum(b==0))
sel = np.logical_and(a != 0 , b !=0 )
vmin = min(a[sel].min(), b[sel].min())
vmax = max(a[sel].max(), b[sel].max())
plot([vmin, vmax], [vmin, vmax], color='tomato')
plot( a[sel], b[sel], '.',ms=1, color='C0')
ax = gca()
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_xlabel("|F| before stage2", fontsize=15)
ax.set_ylabel("|F| after stage2", fontsize=15)
ax.tick_params(labelsize=12)

figure()
res = fobs.d_spacings().data().as_numpy_array()
perc_change = np.abs(b[sel] - a[sel])/a[sel]* 100.
plot(res[sel], perc_change, '.', ms=1)
ax = gca()
ax.set_yscale("log")
ax.invert_xaxis()
ax.tick_params(labelsize=12)
ax.set_xlabel("resolution $\AA$", fontsize=15)
ax.set_ylabel("% change in amplitude", fontsize=15)

#D = np.array(all_data)

show()

