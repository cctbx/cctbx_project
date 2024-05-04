from __future__ import division, print_function
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("mtzin", help="input mtz file", type=str)
parser.add_argument("mtzout", help="output mtz file", type=str)
args = parser.parse_args()

# LIBTBX_SET_DISPATCHER_NAME diffBragg.completeF

import numpy as np
from iotbx.reflection_file_reader import any_reflection_file
from dials.array_family import flex
from scipy.interpolate import interp1d
from cctbx import miller

F = any_reflection_file(args.mtzin).as_miller_arrays()[0]
if not F.is_xray_amplitude_array():
    F = F.as_amplitude_array()
assert F.is_xray_amplitude_array()

print("Bin-ID    Res-range    Completeness    #ASU-indices")
F.show_completeness()
d_max,d_min = F.resolution_range()
print("d_min, d_max (Angstrom): ", d_min, d_max)
mset_full = F.build_miller_set(True, d_min=d_min)
mset_full_d = {h: d for h,d in zip(mset_full.d_spacings().indices(), mset_full.d_spacings().data())}
Fmap = {h:val for h,val in zip(F.indices(), F.data())}
xvals = np.array(F.d_spacings().data())
yvals = np.array(F.data())
fill_vals = yvals[np.argmin(xvals)], yvals[np.argmax(xvals)]
I = interp1d(xvals, yvals, fill_value=fill_vals, bounds_error=False)
data = []
for h in mset_full.indices():
    if h not in Fmap:
        d_h = mset_full_d[h]
        amp = I(d_h)
    else:
        amp = Fmap[h]
    data.append(amp)

complete_amps = flex.double(data)
complete_inds = mset_full.indices()
ma = miller.array(mset_full, complete_amps)
if not ma.is_xray_amplitude_array():
    ma = ma.set_observation_type_xray_amplitude()
ma = ma.as_anomalous_array()
assert ma.anomalous_flag()
ma.as_mtz_dataset(column_root_label="F").mtz_object().write(args.mtzout)
