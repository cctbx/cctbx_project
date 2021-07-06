from __future__ import division
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
parser.add_argument("--curvatures", action='store_true')
parser.add_argument("--readout", type=float, default=0)
args = parser.parse_args()

from dxtbx.model.crystal import Crystal
from cctbx import uctbx
from scitbx.matrix import sqr, rec, col
import numpy as np
from scipy.spatial.transform import Rotation
from simtbx.diffBragg import refine_launcher
from simtbx.nanoBragg.nanoBragg_crystal import NBcrystal
from simtbx.nanoBragg.sim_data import SimData
from simtbx.diffBragg import utils
from dxtbx.model import Experiment
from simtbx.nanoBragg import make_imageset
from cctbx_project.simtbx.diffBragg.phil import phil_scope

ucell = (55, 65, 75, 90, 95, 90)
ucell2 = (55.1, 65.2, 74.9, 90, 94.9, 90)
symbol = "P121"

# generate a random raotation
rotation = Rotation.random(num=1, random_state=100)[0]
Q = rec(rotation.as_quat(), n=(4, 1))
rot_ang, rot_axis = Q.unit_quaternion_as_axis_and_angle()

# generate a small perturbation rotation
np.random.seed(1)
perturb_rot_axis = np.random.random(3)
perturb_rot_axis /= np.linalg.norm(perturb_rot_axis)
perturb_rot_ang = 0.15  # degree random perturbtation

# make the ground truth crystal:
a_real, b_real, c_real = sqr(uctbx.unit_cell(ucell).orthogonalization_matrix()).transpose().as_list_of_lists()
C = Crystal(a_real, b_real, c_real, symbol)
C.rotate_around_origin(rot_axis, rot_ang)

a2_real, b2_real, c2_real = sqr(uctbx.unit_cell(ucell2).orthogonalization_matrix()).transpose().as_list_of_lists()
C2 = Crystal(a2_real, b2_real, c2_real, symbol)
C2.rotate_around_origin(rot_axis, rot_ang)
assert np.allclose(C2.get_U(), C.get_U())
C2.rotate_around_origin(col(perturb_rot_axis), perturb_rot_ang)

# Setup the simulation and create a realistic image
# with background and noise
# <><><><><><><><><><><><><><><><><><><><><><><><><>
nbcryst = NBcrystal()
nbcryst.dxtbx_crystal = C   # simulate ground truth
nbcryst.thick_mm = 0.1
nbcryst.Ncells_abc = 12, 12, 12

SIM = SimData()
SIM.detector = SimData.simple_detector(150, 0.1, (513, 512))
SIM.crystal = nbcryst
SIM.instantiate_diffBragg(oversample=0, auto_set_spotscale=True)
SIM.D.default_F = 0
SIM.D.F000 = 0
SIM.D.progress_meter = False
SIM.water_path_mm = 0.005
SIM.air_path_mm = 0.1
SIM.add_air = True
SIM.add_Water = True
SIM.include_noise = True
SIM.D.add_diffBragg_spots()
spots = SIM.D.raw_pixels.as_numpy_array()
SIM._add_background()
SIM.D.readout_noise_adu=args.readout
SIM._add_noise()
# This is the ground truth image:
img = SIM.D.raw_pixels.as_numpy_array()
SIM.D.raw_pixels *= 0

## Simulate the perturbed image for comparison
#SIM.D.Bmatrix = C2.get_B()
#SIM.D.Umatrix = C2.get_U()
#SIM.D.add_diffBragg_spots()
#SIM._add_background()
#SIM._add_noise()
## Perturbed image:
#img_pet = SIM.D.raw_pixels.as_numpy_array()
#SIM.D.raw_pixels *= 0

E = Experiment()
E.detector = SIM.detector
E.beam = SIM.D.beam
E.crystal = C2  # intentionally set the wrong xtal model
E.imageset = make_imageset([img], E.beam, E.detector)

refls = utils.refls_from_sims([img], E.detector, E.beam, thresh=20)

P = phil_scope.extract()
P.roi.shoebox_size = 20
P.roi.reject_edge_reflections = False
P.refiner.refine_Umatrix = [1]
P.refiner.refine_bg = [0]
P.refiner.refine_Bmatrix = [1]
P.refiner.sensitivity.unitcell = [1, 1, 1, 0.1, 0.1, 0.1]
P.refiner.sensitivity.rotXYZ = [.1, .1, .05]
P.refiner.sensitivity.spot_scale = 0.05
P.refiner.max_calls = [10000]
P.refiner.tradeps = 1e-10
# NOTE RUC.gtol = .9
# NOTE RUC.trad_conv = True  #False
# NOTE RUC.drop_conv_max_eps = 1e-9
P.refiner.curvatures = args.curvatures
P.refiner.use_curvatures_threshold = 0
P.refiner.poissononly = args.readout == 0
P.refiner.verbose = True
P.refiner.big_dump = False
P.refiner.sigma_r = SIM.D.readout_noise_adu
P.refiner.adu_per_photon = SIM.D.quantum_gain
P.refiner.init.ncells_abc = 12, 12, 12
P.refiner.sensitivity.ncells_abc = [1,1,1]
P.simulator.crystal.has_isotropic_ncells = True
P.simulator.crystal.ncells_abc = 12, 12, 12
P.simulator.init_scale = SIM.D.spot_scale
P.simulator.beam.size_mm = SIM.beam.size_mm

#assert RUC.all_ang_off[0] < 0.005
RUC = refine_launcher.local_refiner_from_parameters(refls, E, P, miller_data=SIM.crystal.miller_array)
Ccorrect = RUC.get_corrected_crystal(i_shot=0)
misset = utils.compare_with_ground_truth(*C.get_real_space_vectors(), dxcryst_models=[Ccorrect], symbol=symbol)
assert misset[0] < 0.005
print(misset, "misset with ground truth")
print("OK")
