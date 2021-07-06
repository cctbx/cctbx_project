from __future__ import division
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
parser.add_argument("--curvatures", action='store_true')
args = parser.parse_args()

from dxtbx.model.crystal import Crystal
from cctbx import uctbx
from scitbx.matrix import sqr, rec
from scipy.spatial.transform import Rotation
import numpy as np

from simtbx.nanoBragg.nanoBragg_crystal import NBcrystal
from simtbx.nanoBragg.sim_data import SimData
from simtbx.diffBragg import utils


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
ncells_def_gt = 10,11,12
ncells_abc_gt = 23,25,15
nbcryst.Ncells_abc = ncells_abc_gt
nbcryst.Ncells_def = ncells_def_gt

print("Ground truth ncells DEF abc=%f %f %f" % (ncells_def_gt)) 

# generate the ground truth image
SIM = SimData()
SIM.detector = SimData.simple_detector(200, 0.1, (1024, 1024))
SIM.crystal = nbcryst
SIM.instantiate_diffBragg(oversample=0,auto_set_spotscale=True)
SIM.D.progress_meter = False
SIM.water_path_mm = 0.005
SIM.air_path_mm = 0.1
SIM.add_air = True
SIM.add_Water = True
SIM.include_noise = True
SIM.D.verbose = 2
SIM.D.add_diffBragg_spots()
print("DONE")
spots = SIM.D.raw_pixels.as_numpy_array()
SIM._add_background()
SIM._add_noise()
# This is the ground truth image:
img = SIM.D.raw_pixels.as_numpy_array()
SIM.D.raw_pixels *= 0

# perturb and then
# simulate the perturbed image for comparison
Ncells_def_2 = 7,8,9

SIM.crystal.Ncells_def = Ncells_def_2
SIM.D.Ncells_def = Ncells_def_2
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
P.refiner.refine_ncells_def = [1]
P.refiner.refine_ncells = [1]
P.refiner.ranges.ncells_abc = [3,100,3,100,3,100]
P.refiner.refine_spot_scale = [1]
P.refiner.max_calls = [100]
P.refiner.tradeps = 1e-10
P.refiner.curvatures = args.curvatures
P.refiner.use_curvatures_threshold = 15
P.refiner.poissononly = False
P.refiner.verbose = 2
P.refiner.big_dump = False
P.refiner.sigma_r = SIM.D.readout_noise_adu
P.refiner.adu_per_photon = SIM.D.quantum_gain
Ncells_abc_2 = 10,10,10 # perturb the Ncells abc
P.simulator.crystal.ncells_abc = Ncells_abc_2
P.simulator.crystal.ncells_def = Ncells_def_2
offset_scale = 100
P.simulator.init_scale = SIM.D.spot_scale*offset_scale
P.simulator.beam.size_mm = SIM.beam.size_mm

RUC = refine_launcher.local_refiner_from_parameters(refls, E, P, miller_data=SIM.crystal.miller_array)
Ncells_def_opt = RUC._get_ncells_def_vals(0)
print("Ncells def guess:", Ncells_def_2)
print("Ncells def optimized: ", Ncells_def_opt)
print("Ncells def GT: ", ncells_def_gt)

Ncells_abc_opt = RUC._get_m_val(0)
print("Ncells abc guess:", Ncells_abc_2)
print("Ncells abc optimized: ", Ncells_abc_opt)
print("Ncells abc GT: ", ncells_abc_gt)

deviation = np.sqrt(sum((np.array(Ncells_def_opt) - np.array(ncells_def_gt))**2))
assert deviation < .3, deviation
deviation = np.sqrt(sum((np.array(Ncells_abc_opt) - np.array(ncells_abc_gt))**2))
assert deviation < .3, deviation
assert np.abs(offset_scale - 1/(RUC._get_spot_scale(0)**2))  / offset_scale  < .01

print("OK!")
