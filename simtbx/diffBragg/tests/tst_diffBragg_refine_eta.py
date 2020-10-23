from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action="store_true")
parser.add_argument("--finitediff", action="store_true")
args = parser.parse_args()

from dxtbx.model.crystal import Crystal
from cctbx import uctbx
from scitbx.matrix import rec, col
import numpy as np
from scipy.spatial.transform import Rotation
from scitbx.matrix import sqr
from simtbx.diffBragg.nanoBragg_crystal import nanoBragg_crystal
from simtbx.diffBragg.sim_data import SimData
from simtbx.diffBragg.utils import fcalc_from_pdb
import pylab as plt
import sys


ucell = (79, 79, 38, 90, 90, 90)
symbol = "P43212"
N_MOS_DOMAINS = 100
MOS_SPREAD = 1
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
nbcryst = nanoBragg_crystal()
nbcryst.dxtbx_crystal = C   # simulate ground truth
nbcryst.thick_mm = 0.1
nbcryst.Ncells_abc = Ncells_gt  # ground truth Ncells
nbcryst.mos_spread_deg = MOS_SPREAD
nbcryst.n_mos_domains = N_MOS_DOMAINS
nbcryst.miller_array = miller_array_GT
print("Ground truth ncells = %f" % (nbcryst.Ncells_abc[0]))

# ground truth detector
DET_gt = SimData.simple_detector(150, 0.177, (600, 600))

# initialize the simulator
SIM = SimData()
SIM.Umats_method = 1
SIM.detector = DET_gt
SIM.crystal = nbcryst
SIM.instantiate_diffBragg(oversample=1, verbose=0)
if args.finitediff:
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
SIM.D.add_diffBragg_spots()

if args.finitediff:
  deriv = SIM.D.get_derivative_pixels(eta_diffBragg_id).as_numpy_array()
  SPOTS = SIM.D.raw_pixels_roi.as_numpy_array()
else:
  SPOTS = SIM.D.raw_pixels.as_numpy_array()

SIM.D.readout_noise_adu = 1
SIM._add_background()
SIM._add_noise()

# This is the ground truth image:
img = SIM.D.raw_pixels.as_numpy_array()
SIM.D.raw_pixels_roi *= 0
SIM.D.raw_pixels *= 0

if args.finitediff:
  all_errors = []
  all_shifts = []
  for finite_diff_step in [1, 2, 4, 8, 16]:
    # update Umats to do finite difference test
    delta_eta = 0.0001*finite_diff_step
    SIM.update_umats(MOS_SPREAD + delta_eta, N_MOS_DOMAINS)
    SIM.D.add_diffBragg_spots()

    img_forward = SIM.D.raw_pixels_roi.as_numpy_array()
    SIM.D.raw_pixels_roi *= 0
    SIM.D.raw_pixels *= 0
    fdiff = (img_forward - SPOTS) / delta_eta
    bragg = SPOTS > 1e-2
    error = (np.abs(fdiff[bragg] - deriv[bragg])).mean()
    all_errors.append(error)
    all_shifts.append(delta_eta)

  from scipy.stats import linregress
  l = linregress(all_shifts, all_errors)
  print("finite diff l.rvalue=%10.7g" % l.rvalue)
  assert l.rvalue > .99, "%f" % l.rvalue
  assert l.slope > 0, "%f" % l.slope
  assert l.pvalue < 1e-6, "%f" % l.pvalue
  print("OK")
  if args.plot:
    plt.figure()
    plt.plot(all_shifts, all_errors, 'o')
    plt.show()
  sys.exit()


from simtbx.diffBragg import utils
from dxtbx.model import Experiment
from simtbx.nanoBragg import make_imageset
from cctbx_project.simtbx.diffBragg.phil import phil_scope
from simtbx.diffBragg import refine_launcher

E = Experiment()
E.detector = DET_gt
E.beam = SIM.beam.nanoBragg_constructor_beam
E.crystal = C
E.imageset = make_imageset([img], E.beam, E.detector)

refls = utils.refls_from_sims([SPOTS], E.detector, E.beam, thresh=2000)
print("Found %d spots for refinement" % len(refls))

P = phil_scope.extract()
P.roi.shoebox_size = 25
P.roi.reject_edge_reflections = False
P.refiner.refine_eta = [1]
P.simulator.crystal.mosaicity = 0.1
P.simulator.crystal.num_mosaicity_samples = N_MOS_DOMAINS
P.simulator.beam.size_mm = SIM.beam.size_mm
total_flux = SIM.beam.spectrum[0][1]
print("total flux %f " % total_flux)
P.simulator.total_flux = total_flux
P.simulator.oversample = 1
P.refiner.max_calls = [20]
P.refiner.tradeps = 1e-7
P.refiner.curvatures = False
P.refiner.use_curvatures_threshold = 0
P.refiner.poissononly = False
P.refiner.verbose = True
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
eta_refined = RUC._get_eta(i_shot=0)
print("Refined eta: %10.7f" % eta_refined)
print("GT eta: %10.7f" % MOS_SPREAD)
assert abs(eta_refined - MOS_SPREAD) < 0.001

print("OK")
