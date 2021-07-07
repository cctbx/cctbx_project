from __future__ import division
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action="store_true")
parser.add_argument("--finitediff", action="store_true")
parser.add_argument("--cuda", action="store_true")
args = parser.parse_args()

from dxtbx.model.crystal import Crystal
from cctbx import uctbx
from scitbx.matrix import rec, col
import numpy as np
from scipy.spatial.transform import Rotation
from scitbx.matrix import sqr
from simtbx.nanoBragg.nanoBragg_crystal import NBcrystal
from simtbx.nanoBragg.sim_data import SimData
from simtbx.diffBragg.utils import fcalc_from_pdb
import pylab as plt

SAUSAGE_ID = 20

ucell = (79, 79, 38, 90, 90, 90)
symbol = "P43212"
sausages_diffBragg_id = 20

miller_array_GT = fcalc_from_pdb(resolution=2, wavelength=1, algorithm='fft', ucell=ucell, symbol=symbol)
Ncells_gt = 15, 15, 15

np.random.seed(3142019)
# generate a random rotation
rotation = Rotation.random(num=1, random_state=112)[0]
Q = rec(rotation.as_quat(), n=(4, 1))
rot_ang, rot_axis = Q.unit_quaternion_as_axis_and_angle()

a_real, b_real, c_real = sqr(uctbx.unit_cell(ucell).orthogonalization_matrix()).transpose().as_list_of_lists()
x = col((-1, 0, 0))
y = col((0, -1, 0))
z = col((0, 0, -1))
rx, ry, rz = np.random.uniform(-180, 180, 3)
RX = x.axis_and_angle_as_r3_rotation_matrix(rx, deg=True)
RY = y.axis_and_angle_as_r3_rotation_matrix(ry, deg=True)
RZ = z.axis_and_angle_as_r3_rotation_matrix(rz, deg=True)
M = RX * RY * RZ
a_real = M * col(a_real)
b_real = M * col(b_real)
c_real = M * col(c_real)
C = Crystal(a_real, b_real, c_real, symbol)
C.rotate_around_origin(rot_axis, rot_ang)

# Setup the simulation and create a realistic image
# with background and noise
# <><><><><><><><><><><><><><><><><><><><><><><><><>
nbcryst = NBcrystal(init_defaults=True)
nbcryst.dxtbx_crystal = C   # simulate ground truth
nbcryst.thick_mm = 0.1
nbcryst.Ncells_abc = Ncells_gt  # ground truth Ncells
nbcryst.mos_spread_deg = 0  # MOS_SPREAD
nbcryst.n_mos_domains = 1
nbcryst.miller_array = miller_array_GT
print("Ground truth ncells = %f" % (nbcryst.Ncells_abc[0]))

# sausage 1:
x = col((-1, 0, 0))
y = col((0, -1, 0))
z = col((0, 0, -1))
thetaX, thetaY, thetaZ = 0.1, 0.2, 0.3
RX = x.axis_and_angle_as_r3_rotation_matrix(thetaX, deg=True)
RY = y.axis_and_angle_as_r3_rotation_matrix(thetaY, deg=True)
RZ = z.axis_and_angle_as_r3_rotation_matrix(thetaZ, deg=True)
RXYZ = RX*RY*RZ

# sausage 2
thetaX2, thetaY2, thetaZ2 = -0.3, 0.05, 0.01
RX2 = x.axis_and_angle_as_r3_rotation_matrix(thetaX2, deg=True)
RY2 = y.axis_and_angle_as_r3_rotation_matrix(thetaY2, deg=True)
RZ2 = z.axis_and_angle_as_r3_rotation_matrix(thetaZ2, deg=True)
RXYZ2 = RX2*RY2*RZ2

# ground truth detector
img_sh = 600, 600
DET_gt = SimData.simple_detector(150, 0.177, img_sh)

# initialize the simulator
SIM = SimData()
SIM.detector = DET_gt
SIM.crystal = nbcryst
SIM.instantiate_diffBragg(oversample=1, verbose=0, auto_set_spotscale=True)
#if args.finitediff:
#  SIM.D.refine(eta_diffBragg_id)
#  SIM.D.initialize_managers()
Uorig = sqr(SIM.D.Umatrix)

SIM.D.spot_scale = 100000
SIM.D.default_F = 0
SIM.D.progress_meter = False
SIM.water_path_mm = 0.15
SIM.air_path_mm = 0.1
SIM.add_air = True
SIM.add_water = True
SIM.include_noise = True
SIM.D.use_cuda = args.cuda
SIM.D.printout_pixel_fastslow = 100, 100
SIM.D.add_diffBragg_spots()
SPOTS = SIM.D.raw_pixels_roi.as_numpy_array()
sausages_scale = 0.5, 0.1
sausage_matrices = [RXYZ, RXYZ2]
sausage_spots = []
print(Uorig.elems)
for i in range(2):
    U = sausage_matrices[i]*Uorig
    SIM.D.raw_pixels_roi *= 0
    print(U.elems)
    SIM.D.Umatrix = U.elems
    SIM.D.add_diffBragg_spots()
    sausage_spots.append(sausages_scale[i]*SIM.D.raw_pixels_roi.as_numpy_array())

SIM.D.Umatrix = Uorig
SIM.D.raw_pixels_roi *= 0
NUM_SAUSAGES = 3
SIM.D.update_number_of_sausages(NUM_SAUSAGES)
SIM.D.refine(SAUSAGE_ID)
SIM.D.printout_pixel_fastslow = 100, 100
SIM.D.add_diffBragg_spots()

derivs = SIM.D.get_sausage_derivative_pixels()
if np.allclose( derivs[0], derivs[4]):
  print(1)
else:
  print(2)

SIM.D.raw_pixels_roi *= 0
SIM.D.raw_pixels *= 0

# set the sausages directly in diffBragg
from scitbx.array_family import flex
rotX = flex.double([0, thetaX, thetaX2])*np.pi/180
rotY = flex.double([0, thetaY, thetaY2])*np.pi/180
rotZ = flex.double([0, thetaZ, thetaZ2])*np.pi/180
sausage_scale, sausage_scale2 = sausages_scale
scales = flex.sqrt(flex.double([1, sausage_scale, sausage_scale2]))
SIM.D.set_sausages(rotX, rotY, rotZ, scales)

SIM.D.printout_pixel_fastslow = (100, 100)
SIM.D.add_diffBragg_spots()
sausage_img = SIM.D.raw_pixels_roi.as_numpy_array().reshape(img_sh)

SPOTS2, SPOTS3 = sausage_spots
all_spots = SPOTS+SPOTS2+SPOTS3
manual_adding = all_spots.reshape(img_sh)
assert np.allclose(sausage_img, manual_adding)
print("OK!")

SIM.D.raw_pixels = flex.double(SPOTS+SPOTS2+SPOTS3)
SIM.D.readout_noise_adu = 1
SIM._add_background()
SIM._add_noise()
img = SIM.D.raw_pixels.as_numpy_array()
SIM.D.printout = False

derivs = [d.as_numpy_array() for d in derivs]

if args.finitediff:
  for sausage_idx in [0, 1, 2]:
    for param_idx in range(4):  # there are 4 parameters per sausage
      deriv = derivs[sausage_idx * 4 + param_idx]
      all_errors = []
      all_shifts = []
      for finite_diff_step in [1, 2, 4, 8, 16]:
        delta = 0.0001*finite_diff_step

        rotX = flex.double([0, thetaX, thetaX2]) * np.pi / 180
        rotY = flex.double([0, thetaY, thetaY2]) * np.pi / 180
        rotZ = flex.double([0, thetaZ, thetaZ2]) * np.pi / 180
        scales = flex.sqrt(flex.double([1, sausage_scale, sausage_scale2]))

        if param_idx == 0:
          rotX[sausage_idx] += delta
        elif param_idx == 1:
          rotY[sausage_idx] += delta
        elif param_idx == 2:
          rotZ[sausage_idx] += delta
        else:
          scales[sausage_idx] += delta

        SIM.D.set_sausages(rotX, rotY, rotZ, scales)

        SIM.D.raw_pixels_roi *= 0
        SIM.D.add_diffBragg_spots()

        img_forward = SIM.D.raw_pixels_roi.as_numpy_array()
        fdiff = (img_forward - sausage_img.ravel()) / delta
        bragg = SPOTS > 1e-2

        error = (np.abs(fdiff[bragg] - deriv[bragg]))
        print("MAX ERR: %f" % error.max())
        error = error.mean()
        all_errors.append(error)
        all_shifts.append(delta)
      if args.plot:
        plt.figure()
        plt.plot(all_shifts, all_errors, 'o')
        plt.show()

      from scipy.stats import linregress
      l = linregress(all_shifts, all_errors)
      print("finite diff l.rvalue=%10.7g" % l.rvalue)
      assert l.rvalue > .99, "%f" % l.rvalue
      assert l.slope > 0, "%f" % l.slope
      assert l.pvalue < 1e-6, "%f" % l.pvalue
      print("OK")


from simtbx.diffBragg import utils
from dxtbx.model import Experiment
from simtbx.nanoBragg import make_imageset
from cctbx_project.simtbx.diffBragg.phil import phil_scope
from simtbx.diffBragg import refine_launcher

E = Experiment()
E.detector = DET_gt
E.beam = SIM.beam.nanoBragg_constructor_beam
E.crystal = C
E.imageset = make_imageset([img.reshape(img_sh)], E.beam, E.detector)

refls = utils.refls_from_sims([all_spots.reshape(img_sh)], E.detector, E.beam, thresh=100)
print("Found %d spots for refinement" % len(refls))

P = phil_scope.extract()
P.roi.shoebox_size = 25
P.roi.reject_edge_reflections = False
P.refiner.refine_blueSausages = [1]
num_sausages_guess = 3
P.simulator.crystal.num_sausages = num_sausages_guess
#P.refiner.init.sausages = [0, 0, 0, 1] * num_sausages_guess
P.simulator.beam.size_mm = SIM.beam.size_mm
total_flux = SIM.beam.spectrum[0][1]
print("total flux %f " % total_flux)
P.simulator.total_flux = total_flux
P.simulator.oversample = 1
P.refiner.max_calls = [100]
P.refiner.tradeps = 1e-7
P.refiner.curvatures = False
P.refiner.use_curvatures_threshold = 0
P.refiner.poissononly = False
P.refiner.verbose = 0
P.refiner.big_dump = False
P.refiner.sigma_r = SIM.D.readout_noise_adu
if args.plot:
  P.refiner.plot.display = True
P.refiner.adu_per_photon = SIM.D.quantum_gain
P.simulator.crystal.has_isotropic_ncells = True
P.simulator.init_scale = SIM.D.spot_scale
P.simulator.crystal.ncells_abc = Ncells_gt
P.refiner.ncells_mask = "111"

RUC = refine_launcher.local_refiner_from_parameters(refls, E, P, miller_data=SIM.crystal.miller_array)
sausage_params = RUC._get_sausage_parameters(0)


rotX = flex.double(sausage_params[0::4])
rotY = flex.double(sausage_params[1::4])
rotZ = flex.double(sausage_params[2::4])
scales = flex.double(sausage_params[3::4])
SIM.D.update_number_of_sausages(num_sausages_guess)
SIM.D.set_sausages(rotX, rotY, rotZ, scales)
SIM.D.raw_pixels_roi*= 0
SIM.D.add_diffBragg_spots()
img_after_refinement = SIM.D.raw_pixels_roi.as_numpy_array().reshape(img_sh)

# check pre-refinement image
rotX = flex.double([0]*num_sausages_guess)
rotY = flex.double([0]*num_sausages_guess)
rotZ = flex.double([0]*num_sausages_guess)
scales = flex.double([1]*num_sausages_guess)
SIM.D.set_sausages(rotX, rotY, rotZ, scales)
SIM.D.raw_pixels_roi*= 0
SIM.D.add_diffBragg_spots()
img_before_refinement =SIM.D.raw_pixels_roi.as_numpy_array().reshape(img_sh)

from scipy.stats import pearsonr
a = pearsonr(manual_adding.ravel(), img_after_refinement.ravel())[0]
b = pearsonr(manual_adding.ravel(), img_before_refinement.ravel())[0]
print("Before refinement img/model pearsonR: %f" %b)
print("After refinement img/model pearsonR: %f" %a)
assert a > b, "a/b %f/%f" % (a, b)
assert a > 0.998, "a/b %f/%f" % (a, b)  # this is with 100 iterations - more iterations will give higher correlation (max_calls)

# first test results
# a = 0.998
# b = 0.92

print("OK!")
