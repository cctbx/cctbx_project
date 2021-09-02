from __future__ import division
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action="store_true")
parser.add_argument("--curvatures", action="store_true")
parser.add_argument("--cuda", action="store_true")
parser.add_argument("--aniso", type=int, choices=[0,1,2], default=None)
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


ucell = (79, 79, 38, 90, 90, 90)
symbol = "P43212"
N_MOS_DOMAINS = 100
MOS_SPREAD = 1
ANISO_MOS_SPREAD = 0.5, 0.75, 1
eta_diffBragg_id = 19

miller_array_GT = fcalc_from_pdb(resolution=2, wavelength=1, algorithm='fft', ucell=ucell, symbol=symbol)
Ncells_gt = 15, 15, 15

np.random.seed(3142019)
# generate a random rotation
rotation = Rotation.random(num=1, random_state=100)[0]
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
nbcryst.mos_spread_deg = MOS_SPREAD
if args.aniso is not None:
  nbcryst.anisotropic_mos_spread_deg = ANISO_MOS_SPREAD
  assert nbcryst.has_anisotropic_mosaicity
else:
  assert not nbcryst.has_anisotropic_mosaicity

nbcryst.n_mos_domains = N_MOS_DOMAINS
nbcryst.miller_array = miller_array_GT
print("Ground truth ncells = %f" % (nbcryst.Ncells_abc[0]))

# ground truth detector
DET_gt = SimData.simple_detector(150, 0.177, (600, 600))

# initialize the simulator
SIM = SimData()
if args.aniso is None:
  SIM.Umats_method = 2
else:
  SIM.Umats_method = 3
SIM.detector = DET_gt
SIM.crystal = nbcryst
SIM.instantiate_diffBragg(oversample=1, verbose=0)
SIM.D.refine(eta_diffBragg_id)
SIM.D.initialize_managers()
SIM.D.spot_scale = 100000
SIM.D.default_F = 0
SIM.D.progress_meter = False
SIM.water_path_mm = 0.15
SIM.air_path_mm = 0.1
SIM.add_air = True
SIM.add_water = True
SIM.include_noise = True
SIM.D.use_cuda = args.cuda
SIM.D.compute_curvatures = args.curvatures
SIM.D.add_diffBragg_spots()

if args.aniso is None:
  deriv = SIM.D.get_derivative_pixels(eta_diffBragg_id).as_numpy_array()
else:
  deriv = SIM.D.get_aniso_eta_deriv_pixels()[args.aniso].as_numpy_array()

if args.curvatures:
  if args.aniso is None:
    deriv2 = SIM.D.get_second_derivative_pixels(eta_diffBragg_id).as_numpy_array()
  else:
    deriv2 = SIM.D.get_aniso_eta_second_deriv_pixels()[args.aniso].as_numpy_array()
SPOTS = SIM.D.raw_pixels_roi.as_numpy_array()

SIM.D.readout_noise_adu = 1
SIM._add_background()
SIM._add_noise()

# This is the ground truth image:
img = SIM.D.raw_pixels.as_numpy_array()
SIM.D.raw_pixels_roi *= 0
SIM.D.raw_pixels *= 0

all_errors = []
all_shifts = []
all_errors2 = []
all_shifts2 = []
for finite_diff_step in [1, 2, 4, 8, 16]:
  # update Umats to do finite difference test
  delta_eta = 0.001*finite_diff_step

  if args.aniso is not None:
    eta_update = list(ANISO_MOS_SPREAD)
    eta_update[args.aniso] = eta_update[args.aniso]+ delta_eta
    crystal = nbcryst.dxtbx_crystal
  else:
    eta_update = MOS_SPREAD + delta_eta
    crystal = None

  SIM.update_umats(eta_update, N_MOS_DOMAINS, crystal=crystal)
  SIM.D.add_diffBragg_spots()

  img_forward = SIM.D.raw_pixels_roi.as_numpy_array()
  SIM.D.raw_pixels_roi *= 0
  SIM.D.raw_pixels *= 0
  fdiff = (img_forward - SPOTS) / delta_eta
  bragg = SPOTS > 1e-2
  error = (np.abs(fdiff[bragg] - deriv[bragg])).mean()
  all_errors.append(error)
  all_shifts.append(delta_eta)
  if args.curvatures:

    if args.aniso is not None:
      eta_update = list(ANISO_MOS_SPREAD)
      eta_update[args.aniso] = eta_update[args.aniso] - delta_eta
      crystal = nbcryst.dxtbx_crystal
    else:
      eta_update = MOS_SPREAD - delta_eta
      crystal= None

    SIM.update_umats(eta_update, N_MOS_DOMAINS, crystal=crystal)

    all_shifts2.append(delta_eta ** 2)

    SIM.D.raw_pixels_roi *= 0
    SIM.D.raw_pixels *= 0
    SIM.D.add_diffBragg_spots()
    img_backward = SIM.D.raw_pixels_roi.as_numpy_array()

    fdiff2 = (img_forward - 2*SPOTS + img_backward) / delta_eta/ delta_eta

    second_deriv = deriv2
    error2 = (np.abs(fdiff2[bragg] - second_deriv[bragg]) / 1).mean()
    all_errors2.append(error2)

  print("\n\n<><><><><><><><>\n\tError:", error, "shift:", delta_eta)
  if args.curvatures:
    print("\terror2=%f, step=%f\n<><><><><><><><>\n\n" % (error2, delta_eta))

from scipy.stats import linregress
l = linregress(all_shifts, all_errors)
print("finite diff l.rvalue=%10.7g" % l.rvalue)
if args.plot:
  plt.figure()
  plt.plot(all_shifts, all_errors, 'o')
  plt.show()
  if args.curvatures:
    plt.plot(all_shifts2, all_errors2, 'o')
    plt.title("second finite diff error")
    plt.xlabel("delta eta")
    plt.show()

assert l.rvalue > .99, "%f" % l.rvalue
assert l.slope > 0, "%f" % l.slope
assert l.pvalue < 1e-6, "%f" % l.pvalue
if args.curvatures:
  l = linregress(all_shifts2, all_errors2)
  assert l.rvalue > .9999  # this is definitely a line!
  assert l.slope > 0
  assert l.pvalue < 1e-6
print("OK")
