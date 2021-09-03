from __future__ import division
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
parser.add_argument("--curvatures", action='store_true')
parser.add_argument("--readout", type=float, default=0)
parser.add_argument("--perturb", choices=["G", "Nabc", "detz_shift", "crystal"], type=str, nargs="+", default=["crystal"] )
args = parser.parse_args()

from dxtbx.model.crystal import Crystal
from cctbx import uctbx
from scitbx.matrix import sqr, rec, col
import numpy as np
from scipy.spatial.transform import Rotation
from simtbx.nanoBragg.nanoBragg_crystal import NBcrystal
from simtbx.nanoBragg.sim_data import SimData
from simtbx.diffBragg import utils
from dxtbx.model import Experiment
from simtbx.nanoBragg import make_imageset
from simtbx.diffBragg.phil import hopper_phil, philz
from libtbx.phil import parse

phil_scope = parse(hopper_phil+philz)

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
nbcryst.Ncells_abc = 12, 12, 11
nbcryst.isotropic_ncells = False

SIM = SimData(use_default_crystal=True)
#SIM.detector = SimData.simple_detector(150, 0.1, (513, 512))
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

P = phil_scope.extract()
E = Experiment()

if "G" in args.perturb:
    P.init.G = SIM.D.spot_scale*10
else:
    P.init.G = SIM.D.spot_scale

if "crystal" in args.perturb:
    E.crystal = C2
else:
    E.crystal = C

if "Nabc" in args.perturb:
    P.init.Nabc = 20,20,20
else:
    P.init.Nabc = SIM.crystal.Ncells_abc

if "detz_shift" in args.perturb:
    P.init.detz_shift = 1
else:
    P.init.detz_shift = 0

E.detector = SIM.detector
E.beam = SIM.D.beam
E.imageset = make_imageset([img], E.beam, E.detector)
refls = utils.refls_from_sims([img], E.detector, E.beam, thresh=18)
print("%d REFLS" % len(refls))
utils.refls_to_q(refls, E.detector, E.beam, update_table=True)
utils.refls_to_hkl(refls, E.detector, E.beam, E.crystal, update_table=True)

P.roi.shoebox_size = 20
P.relative_tilt = False
P.roi.fit_tilt = False
P.roi.pad_shoebox_for_background_estimation=10
P.roi.reject_edge_reflections = False
P.refiner.sigma_r = SIM.D.readout_noise_adu
P.refiner.adu_per_photon = SIM.D.quantum_gain
P.simulator.init_scale = 1 #SIM.D.spot_scale
P.simulator.beam.size_mm = SIM.beam.size_mm
P.simulator.total_flux = SIM.D.flux
P.use_restraints = False
name = "hopper_refine_%s.mtz" % "-".join(args.perturb) # TODO interface for passing this directly to hopper_utils.refine
SIM.crystal.miller_array.as_mtz_dataset(column_root_label="F").mtz_object().write(name)
P.simulator.structure_factors.mtz_name = name
P.simulator.structure_factors.mtz_column = "F(+),F(-)"
P.niter = 0
P.niter_per_J = 1
P.method="L-BFGS-B"
P.ftol = 1e-10
#P.method="Nelder-Mead"
#P.fix.G = True
#P.fix.Nabc =True
#P.fix.detz_shift=True

import logging
import sys
h = logging.StreamHandler(sys.stdout)
logging.basicConfig(level=logging.DEBUG, handlers=[h])

from simtbx.diffBragg import hopper_utils
Eopt,_, Mod, x = hopper_utils.refine(E, refls, P, return_modeler=True)

G, rotX,rotY, rotZ, Na,Nb,Nc,a,b,c,al,be,ga,detz_shift = hopper_utils.get_param_from_x(x, Mod.SIM)
print("Na, Nb, Nc= %f %f %f" % (Na, Nb, Nc))

# check crystal
Copt = Eopt.crystal
misset, misset_init = utils.compare_with_ground_truth(*C.get_real_space_vectors(), dxcryst_models=[Copt, C2], symbol=symbol)
print(misset_init, "init misset with ground truth")
print(misset, "misset with ground truth")
if "detz_shift" in args.perturb:
    assert misset < 0.0065
else:
    assert misset < 0.005

# check mosaic domain
assert all (np.subtract(nbcryst.Ncells_abc, [Na,Nb,Nc]) < 0.2)

# check spot scale
perc_diff_G = abs(SIM.D.spot_scale - G)/ SIM.D.spot_scale * 100
print("spot scale gt: %f; spot scale opt: %f; percent diff: %f %%" % (SIM.D.spot_scale, G, perc_diff_G))
assert perc_diff_G < 1

# check detz
print("detdist shift %f (should be 0)" % detz_shift)
assert detz_shift < 0.2


ucell_diff_init = np.abs(np.subtract(ucell , ucell2))
ucell_diff = np.abs(np.subtract(ucell , Copt.get_unit_cell().parameters()))


init_dev, init_dev_ang = ucell_diff_init[:3].sum(), ucell_diff_init[-3:].sum()
dev, dev_ang = ucell_diff[:3].sum(), ucell_diff[-3:].sum()
print("initial ucell dev: %f Angstrom; %f degree" % (init_dev, init_dev_ang))
print("optimized ucell dev: %f Angstrom; %f degree" % (dev, dev_ang))
assert dev_ang < init_dev_ang and dev_ang < 0.01
if "detz_shift" not in args.perturb:
    assert dev < init_dev and dev < 0.01

print("OK")
