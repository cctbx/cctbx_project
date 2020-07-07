"""
This test checks the setter and getter for Ncells parameter
"""

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
parser.add_argument("--curvatures", action='store_true')
parser.add_argument("--idx", default=0, help="0=Na, 1=Nb, 2=Nc", type=int, choices=[0,1,2])
args = parser.parse_args()

import numpy as np
import pylab as plt
from scipy.stats import linregress
from scipy.spatial.transform import Rotation
from simtbx.diffBragg import sim_data
from scitbx.matrix import sqr, rec
from dxtbx.model import Detector,Panel
from IPython import embed
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
panel_rot_id = 14  # internal diffBragg index for the panel rotation parameter

S = sim_data.SimData()
S.crystal.dxtbx_crystal = C
S.detector = sim_data.SimData.simple_detector(180, 0.1, (1024, 1024))

DET = S.detector
BEAM = S.beam.nanoBragg_constructor_beam
panel_dict = DET[0].to_dict()
fdet_orig = panel_dict["fast_axis"]
sdet_orig = panel_dict["slow_axis"]
panel = Panel.from_dict(panel_dict)
Dnew = Detector()
Dnew.add_panel(panel)


S.detector = Dnew

S.instantiate_diffBragg(verbose=0, oversample=0)
S.D.spot_scale = 100000
S.D.compute_curvatures = False
S.D.nopolar=False

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
    angle_deg = angle_deg*0.01

    angle_rad = angle_deg*np.pi/180.
    shifts.append(angle_rad)
    S.D.update_dxtbx_geoms(DET, BEAM, 0, angle_rad)
    S.D.raw_pixels *= 0
    S.D.region_of_interest = ((0, 1023), (0, 1023))
    #S.D.printout_pixel_fastslow = 10, 10
    S.D.add_diffBragg_spots()
    img2 = S.D.raw_pixels.as_numpy_array()

    fdiff = (img2 - img) /angle_rad
    error = np.abs(fdiff[bragg] - deriv[bragg]).mean()
    all_error.append(error)
    embed()

    print("error=%f, step=%f radians" % (error, angle_rad))
    #test_ang = np.arccos(np.dot(S.D.fdet_vector, fdet_orig))
    #print("ANGLE between Fdet and starting Fdet = %f, GT ANGLE=%f" % (test_ang*180/np.pi, angle_deg))
    #test_ang = np.arccos(np.dot(S.D.sdet_vector, sdet_orig))
    #print("ANGLE between Sdet and starting Sdet = %f, GT ANGLE=%f" % (test_ang*180/np.pi, angle_deg))
    #print("")

    #if args.plot:
    #    plt.subplot(121)
    #    plt.imshow(fdiff)
    #    plt.title("finite diff")
    #    plt.subplot(122)
    #    plt.imshow(deriv)
    #    plt.title("analytical")
    #    plt.draw()
    #    plt.suptitle("Shift %d / %d"
    #                 % (i_shift + 1, len(perc)))
    #    plt.pause(0.8)

if args.plot:
    #plt.close()
    plt.plot(shifts, all_error, 'o')
    plt.show()
    #if args.curvatures:
    #    plt.plot(shifts2, all_error2, 'o')
    #    plt.show()

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
