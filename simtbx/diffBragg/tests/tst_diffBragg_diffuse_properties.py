"""
This test checks the setter and getter for Ncells parameter
"""
from __future__ import division
from simtbx.kokkos import gpu_instance
kokkos_run = gpu_instance(deviceId = 0)

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--cuda", action="store_true")
parser.add_argument("--plot", action='store_true')
parser.add_argument("--onlyDiffuse", action='store_true')
parser.add_argument("--sigma", default=[1,1,1], type=float, nargs=3)
parser.add_argument("--gamma", default=[100,100,100], type=float, nargs=3)
parser.add_argument("--idx", type=int, default=0, choices=[0,1,2], help="diffuse parameter index (0,1,2 ->a,b,c)")
parser.add_argument("--grad", choices=['sigma','gamma'], default='gamma')

import pylab as plt
args = parser.parse_args()
if args.cuda:
    import os
    os.environ["DIFFBRAGG_USE_CUDA"] = "1"

from simtbx.nanoBragg import sim_data

S = sim_data.SimData(use_default_crystal=True)
det_shape = (1024,1024)
S.detector = sim_data.SimData.simple_detector(180, 0.1, det_shape)
S.instantiate_diffBragg(verbose=0, oversample=0, auto_set_spotscale=True)
#S.D.record_time = True
S.D.spot_scale = 100000
S.D.use_diffuse = True
S.D.diffuse_gamma = tuple(args.gamma)
S.D.diffuse_sigma = tuple(args.sigma)

default_gamma = args.gamma
default_sigma = args.sigma

import numpy as np
assert np.allclose(S.D.diffuse_gamma, args.gamma)
assert np.allclose(S.D.diffuse_sigma, args.sigma)

diffuse_id = 23
S.D.refine(diffuse_id)
S.D.add_diffBragg_spots()
img = S.D.raw_pixels_roi.as_numpy_array()
bragg = img > 1e-1  # select bragg scattering regions

if args.grad == 'gamma':
    derivs_abc = list(map(lambda x: x.as_numpy_array(), S.D.get_diffuse_gamma_derivative_pixels()))
    diff_params = args.gamma
else:
    derivs_abc = list(map(lambda x: x.as_numpy_array(), S.D.get_diffuse_sigma_derivative_pixels()))
    diff_params = args.sigma

S.D.fix(diffuse_id)

from scipy.stats import linregress
perc = 0.001, 0.01, 0.1, 1, 10
for ii in range(3):
    if ii != args.idx:
        continue

    all_error = []
    all_error2 = []
    shifts = []
    shifts2 = []
    gi = diff_params[ii]
    for i_shift, p in enumerate(perc):
        delta_gi = gi*p*0.01
        print(i_shift, delta_gi, gi)
        diff_vals = diff_params*1
        diff_vals[ii] = gi+delta_gi

        if args.grad == 'gamma':
            S.D.diffuse_gamma = tuple(diff_vals)
        else:
            S.D.diffuse_sigma = tuple(diff_vals)

        shifts.append(delta_gi)

        S.D.raw_pixels_roi *= 0
        #S.D.printout_pixel_fastslow = 10, 10
        S.D.add_diffBragg_spots()
        img2 = S.D.raw_pixels_roi.as_numpy_array()

        fdiff = (img2 - img) / delta_gi
        deriv = derivs_abc[ii]
        error = np.abs(fdiff[bragg] - deriv[bragg]).mean()
        all_error.append(error)

        print ("error=%f, step=%f" % (error, delta_gi))

    if args.plot:
        plt.plot(shifts, all_error, 'o')
        plt.show()

    l = linregress(shifts, all_error)
    assert l.rvalue > .9999  # this is definitely a line!
    assert l.slope > 0
    assert l.pvalue < 1e-6

det_sh = 1024, 1024
print("OK")
