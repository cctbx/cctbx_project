"""
This test checks the setter and getter for Ncells parameter
"""
from __future__ import division

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
parser.add_argument("--idx", choices=["odet","sdet", "fdet"], type=str, default="odet")
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
if args.idx == "odet":
    panel_rot_id = 14  # internal diffBragg index for the panel rotation parameter
elif args.idx == "fdet":
    panel_rot_id = 17
else:
    panel_rot_id = 18

S = sim_data.SimData(use_default_crystal=True)
S.crystal.dxtbx_crystal = C
S.detector = sim_data.SimData.simple_detector(180, 0.1, (1024, 1024))

DET = S.detector
BEAM = S.beam.nanoBragg_constructor_beam

S.detector = DET #Dnew

S.instantiate_diffBragg(verbose=0, oversample=0, auto_set_spotscale=True)
S.D.spot_scale = 100000
S.D.compute_curvatures = False
S.D.nopolar = True

S.D.refine(panel_rot_id)
S.D.initialize_managers()
S.D.region_of_interest = ((0, 1023), (0, 1023))
#S.D.printout_pixel_fastslow = 10, 10
S.D.add_diffBragg_spots()
img = S.D.raw_pixels.as_numpy_array()
deriv = S.D.get_derivative_pixels(panel_rot_id).as_numpy_array()
bragg = img > 1e-1  # select bragg scattering regions

# angles in degrees
angles = 0.001, 0.002, 0.004, 0.008, 0.016, 0.032

all_error = []
shifts = []
for i_shift, angle_deg in enumerate(angles):
    scale = 0.01
    if args.idx in ["fdet", "sdet"]:
        scale = 10
    angle_deg = angle_deg*scale

    angle_rad = angle_deg*np.pi/180.
    shifts.append(angle_rad)
    if args.idx == "odet":
        ang1, ang2, ang3 = angle_rad, 0, 0
    elif args.idx == "fdet":
        ang1, ang2, ang3 = 0, angle_rad, 0
    else:
        ang1, ang2, ang3 = 0, 0, angle_rad
    S.D.update_dxtbx_geoms(DET, BEAM, 0, ang1, ang2, ang3)

    S.D.raw_pixels *= 0
    S.D.region_of_interest = ((0, 1023), (0, 1023))
    #S.D.printout_pixel_fastslow = 10, 10
    S.D.add_diffBragg_spots()
    img2 = S.D.raw_pixels.as_numpy_array()

    fdiff = (img2 - img) /angle_rad
    error = np.abs(fdiff[bragg] - deriv[bragg]).mean()
    all_error.append(error)

    print("error=%f, step=%f degrees" % (error, angle_deg))

if args.plot:
    plt.plot(shifts, all_error, 'o')
    plt.show()

l = linregress(shifts, all_error)
assert l.rvalue > .9999  # this is definitely a line!
assert l.slope > 0
assert l.pvalue < 1e-6
#if args.curvatures:
#    l = linregress(shifts2, all_error2)
#    assert l.rvalue > .9999  # this is definitely a line!
#    assert l.slope > 0
#    assert l.pvalue < 1e-6

print("OK!")
