from __future__ import division
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
parser.add_argument("--curvatures", action='store_true')
parser.add_argument("--constrain", action="store_true")
parser.add_argument("--testrangeLower", action="store_true")
parser.add_argument("--testrangeUpper", action="store_true")
args = parser.parse_args()

from dxtbx.model.crystal import Crystal
from cctbx import uctbx
from scitbx.matrix import sqr, rec
from scipy.spatial.transform import Rotation
import numpy as np

from simtbx.nanoBragg.nanoBragg_crystal import NBcrystal
from simtbx.nanoBragg.sim_data import SimData
from simtbx.diffBragg import utils


ucell = (85.2, 96, 124, 90, 105, 90)
symbol = "P121"
ucell = (96, 96, 124, 90, 90, 90)
symbol = "P43212"

# generate a random raotation
rotation = Rotation.random(num=1, random_state=1107)[0]
Q = rec(rotation.as_quat(), n=(4, 1))
rot_ang, rot_axis = Q.unit_quaternion_as_axis_and_angle()

# make the ground truth crystal:
a_real, b_real, c_real = sqr(uctbx.unit_cell(ucell).orthogonalization_matrix()).transpose().as_list_of_lists()
C = Crystal(a_real, b_real, c_real, symbol)
C.rotate_around_origin(rot_axis, rot_ang)

# Setup the simulation and create a realistic image
# with background and noise
# <><><><><><><><><><><><><><><><><><><><><><><><><>
nbcryst = NBcrystal()
nbcryst.dxtbx_crystal = C   # simulate ground truth
nbcryst.thick_mm = 0.1
nbcryst.isotropic_ncells = False
ncells_gt = 19, 23, 15
ncells_gt = 23, 23, 15

nbcryst.Ncells_abc = ncells_gt

print("Ground truth ncells abc=%f %f %f" % (ncells_gt)) #nbcryst.Ncells_abc[0]))

# generate the ground truth image
SIM = SimData()
SIM.detector = SimData.simple_detector(200, 0.1, (1024, 1024))
SIM.crystal = nbcryst
SIM.instantiate_diffBragg(oversample=0, auto_set_spotscale=True)
SIM.D.progress_meter = False
SIM.water_path_mm = 0.005
SIM.air_path_mm = 0.1
SIM.add_air = True
SIM.add_Water = True
SIM.include_noise = True
SIM.D.add_diffBragg_spots()
spots = SIM.D.raw_pixels.as_numpy_array()
SIM._add_background()
SIM._add_noise()
# This is the ground truth image:
img = SIM.D.raw_pixels.as_numpy_array()
SIM.D.raw_pixels *= 0

# perturb and then
# simulate the perturbed image for comparison
Ncells_abc_2 = 26, 26, 26
if args.testrangeLower:
  Ncells_abc_2 = 10, 10, 10

SIM.crystal.Ncells_abc = Ncells_abc_2
SIM.D.set_value(9, Ncells_abc_2[0])
SIM.D.add_diffBragg_spots()
SIM._add_background()
SIM._add_noise()
# Perturbed image:
img_pet = SIM.D.raw_pixels.as_numpy_array()
SIM.D.raw_pixels *= 0

from dxtbx.model import Experiment
from simtbx.nanoBragg import make_imageset
from cctbx_project.simtbx.diffBragg.phil import phil_scope
from simtbx.diffBragg import refine_launcher
E = Experiment()
E.detector = SIM.detector
E.beam = SIM.D.beam
E.crystal = C
E.imageset = make_imageset([img], E.beam, E.detector)

refls = utils.refls_from_sims([spots], E.detector, E.beam, thresh=20)

P = phil_scope.extract()
P.roi.shoebox_size = 20
P.roi.reject_edge_reflections = False
P.refiner.refine_ncells = [1]
P.refiner.max_calls = [10]
P.refiner.tradeps = 1e-10
P.refiner.curvatures = False
P.refiner.use_curvatures_threshold = 0
P.refiner.poissononly = False
P.refiner.verbose = True
P.refiner.big_dump = False
P.refiner.sigma_r = SIM.D.readout_noise_adu
P.refiner.adu_per_photon = SIM.D.quantum_gain
P.simulator.crystal.ncells_abc = Ncells_abc_2  # initial guess
P.simulator.init_scale = SIM.D.spot_scale
P.simulator.beam.size_mm = SIM.beam.size_mm
if args.testrangeLower:
  P.refiner.ranges.ncells_abc = 3,22, 3, 22, 3, 14
elif args.testrangeUpper:
  P.refiner.ranges.ncells_abc = 24,100, 24, 100, 16, 100
if args.constrain:
  P.refiner.ncells_mask = "110"
else:
  P.refiner.ncells_mask = "000"

RUC = refine_launcher.local_refiner_from_parameters(refls, E, P, miller_data=SIM.crystal.miller_array)
Ncells_opt = RUC._get_m_val(0)
print("Ncells guess:", Ncells_abc_2)
print("Ncells optimized: ", Ncells_opt)
print("Ncells GT: ", ncells_gt)
if args.testrangeLower:
  print("... testing bounded refinement ...")
  assert round(Ncells_opt[0]) == round(Ncells_opt[1]) == 22
  assert round(Ncells_opt[2]) == 14
elif args.testrangeUpper:
  print("... testing bounded refinement ...")
  assert round(Ncells_opt[0]) == round(Ncells_opt[1]) == 24
  assert round(Ncells_opt[2]) == 16
else:
  deviation = np.sqrt(sum((np.array(Ncells_opt) - np.array(ncells_gt))**2))
  assert deviation < .3

print("OK!")
