"""
This test checks the setter and getter for Ncells parameter
"""
from __future__ import division

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
parser.add_argument("--curvatures", action='store_true')
parser.add_argument("--idx", default=0, help="0=Nd, 1=Ne, 2=Nf", type=int, choices=[0,1,2])
args = parser.parse_args()

import numpy as np
import pylab as plt
from scipy.stats import linregress
from scipy.spatial.transform import Rotation
from simtbx.nanoBragg import sim_data
from scitbx.matrix import sqr, rec
from cctbx import uctbx
from dxtbx.model import Crystal

ucell = (70, 60, 50, 90.0, 110, 90.0)
symbol = "C121"

a_real, b_real, c_real = sqr(uctbx.unit_cell(ucell).orthogonalization_matrix()).transpose().as_list_of_lists()
C = Crystal(a_real, b_real, c_real, symbol)

# random raotation
rotation = Rotation.random(num=1, random_state=101)[0]
Q = rec(rotation.as_quat(), n=(4, 1))
rot_ang, rot_axis = Q.unit_quaternion_as_axis_and_angle()
C.rotate_around_origin(rot_axis, rot_ang)

S = sim_data.SimData(use_default_crystal=True)
S.crystal.dxtbx_crystal = C
S.detector = sim_data.SimData.simple_detector(180, 0.1, (1024, 1024))
S.instantiate_diffBragg(verbose=0, oversample=0, auto_set_spotscale=True)
S.D.spot_scale = 100000
#S.D.oversample = 5

Ncells_def_id = 21

S.D.Ncells_abc_aniso = 12, 13, 14
S.D.Ncells_def = 10,11,12
Nd, Ne, Nf = S.D.Ncells_def

assert (Nd,Ne,Nf) == (10,11,12)

print("OK")
S.D.isotropic_ncells = False
S.D.refine(Ncells_def_id)
S.D.initialize_managers()
ncells_def_GT = 10,11,12

S.D.Ncells_abc_aniso = (20, 15, 17)
S.D.Ncells_def = ncells_def_GT
S.D.add_diffBragg_spots()

img = S.D.raw_pixels_roi.as_numpy_array()
derivs_abc = S.D.get_ncells_def_derivative_pixels()
derivs_abc = [d.as_numpy_array() for d in derivs_abc]
if args.curvatures:
    derivs2_abc = S.D.get_ncells_def_second_derivative_pixels()
    derivs2_abc = [d.as_numpy_array() for d in derivs2_abc]

perc = 0.001, 0.01, 0.1, 1, 10
bragg = img > 1e-1  # select bragg scattering regions

Nd, Ne, Nf = ncells_def_GT
for i_nc in range(3):
    if i_nc != args.idx:
        continue

    N = ncells_def_GT[i_nc]
    print("Testing finite diff %d / 3, Ncells_def=%d" % (i_nc +1, N))

    all_error = []
    all_error2 = []
    shifts = []
    shifts2 = []
    for i_shift, p in enumerate(perc):
        delta_N = N*p*0.01

        Ncells_def_vals = [Nd, Ne, Nf]
        Ncells_def_vals[i_nc] = N+delta_N

        S.D.Ncells_def = tuple(Ncells_def_vals)
        shifts.append(delta_N)

        S.D.raw_pixels_roi *= 0
        #S.D.printout_pixel_fastslow = 10, 10
        S.D.add_diffBragg_spots()
        img2 = S.D.raw_pixels_roi.as_numpy_array()

        fdiff = (img2 - img) / delta_N
        deriv = derivs_abc[i_nc]
        error = np.abs(fdiff[bragg] - deriv[bragg]).mean()
        all_error.append(error)

        if args.curvatures:
            Ncells_vals = [Nd, Ne, Nf]
            Ncells_vals[i_nc] = N - delta_N
            S.D.Ncells_def = tuple(Ncells_vals)
            shifts2.append(delta_N**2)

            S.D.raw_pixels_roi *= 0
            S.D.add_diffBragg_spots()
            img_backward = S.D.raw_pixels_roi.as_numpy_array()

            fdiff2 = (img2 - 2*img + img_backward) / delta_N / delta_N

            second_deriv = derivs2_abc[i_nc]

            error2 = (np.abs(fdiff2[bragg] - second_deriv[bragg]) / 1).mean()
            all_error2.append(error2)

        print ("error=%f, step=%f" % (error, delta_N))

    if args.plot:
        plt.plot(shifts, all_error, 'o')
        plt.show()
        if args.curvatures:
            plt.plot(shifts2, all_error2, 'o')
            plt.show()

    l = linregress(shifts, all_error)
    assert l.rvalue > .999  # this is definitely a line!
    assert l.slope > 0
    assert l.pvalue < 1e-5
    if args.curvatures:
        l = linregress(shifts2, all_error2)
        assert l.rvalue > .999  # this is definitely a line!
        assert l.slope > 0
        assert l.pvalue < 1e-5

    print("OK!")
