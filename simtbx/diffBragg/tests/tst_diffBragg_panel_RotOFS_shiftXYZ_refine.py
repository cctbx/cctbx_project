from __future__ import division
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
parser.add_argument("--oversample", type=int, default=0)
parser.add_argument("--nopolar", action="store_true")
args = parser.parse_args()

from dxtbx.model.crystal import Crystal
from cctbx import uctbx
from scitbx.matrix import sqr, rec
import numpy as np
from scipy.spatial.transform import Rotation

from simtbx.nanoBragg.nanoBragg_crystal import NBcrystal
from simtbx.nanoBragg.sim_data import SimData
from simtbx.diffBragg import utils

ucell = (85.2, 96, 124, 90, 105, 90)
symbol = "P121"

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
nbcryst.Ncells_abc = 12, 12, 12

SIM = SimData()
SIM.detector = SimData.simple_detector(160, 0.1, (1024, 1024))

# grab the detector node (ground truth)
node = SIM.detector[0]
node_d = node.to_dict()
Origin = node_d["origin"][0], node_d["origin"][1], node_d["origin"][2]
distance = Origin[2]
gt_distance = distance
print("Ground truth originZ=%f" % (SIM.detector[0].get_origin()[2]))


# copy the detector and update the origin
#det2 = deepcopy(SIM.detector)
# alter the detector distance by 2 mm
#node_d["origin"] = Origin[0], Origin[1], Origin[2]+3
#embed()
#det2[0] = Panel.from_dict(node_d)
#print ("Modified originZ=%f" % (det2[0].get_origin()[2]))

SIM.crystal = nbcryst
SIM.instantiate_diffBragg(oversample=args.oversample, auto_set_spotscale=True)
SIM.D.progress_meter = False
SIM.D.verbose = 2 #0 #1
SIM.D.nopolar = args.nopolar
SIM.water_path_mm = 0.005
SIM.air_path_mm = 0.1
SIM.add_air = True
SIM.add_Water = True
SIM.include_noise = True
#SIM.D.spot_scale = 1e8

SIM.D.add_diffBragg_spots()
spots = SIM.D.raw_pixels.as_numpy_array()
SIM._add_background()
SIM._add_noise()
# This is the ground truth image:
img = SIM.D.raw_pixels.as_numpy_array()
SIM.D.raw_pixels *= 0


# get the fast,slow, origin axis of the panel
from scitbx.matrix import col
panel = SIM.detector[0]
Ftru = col(panel.get_fast_axis())
Stru = col(panel.get_slow_axis())
Otru = col(panel.get_origin())

# Simulate the perturbed image for comparison
panel_rot_angO = 0.05
panel_rot_angF = 1
panel_rot_angS = 1
panel_rot_ang_radO = panel_rot_angO * np.pi / 180
panel_rot_ang_radF = panel_rot_angF * np.pi / 180
panel_rot_ang_radS = panel_rot_angS * np.pi / 180
Xshift_mm = 0.05
Yshift_mm = 0.04
Zshift_mm = 0.1

#panel_rot_ang_radS = panel_rot_ang_radF = panel_rot_ang_radO = 0
SIM.D.reference_origin = SIM.detector[0].get_origin()
SIM.D.update_dxtbx_geoms(SIM.detector, SIM.beam.nanoBragg_constructor_beam, 0,
                         panel_rot_ang_radO, panel_rot_ang_radF, panel_rot_ang_radS,
                         Xshift_mm/1000., Yshift_mm/1000., Zshift_mm/1000.)

panel_dict = SIM.detector[0].to_dict()
panel_dict["fast_axis"] = SIM.D.fdet_vector
panel_dict["slow_axis"] = SIM.D.sdet_vector
panel_dict["origin"] = SIM.D.get_origin()

from dxtbx.model import Detector, Panel
perturb_detector = Detector()
perturb_detector.add_panel(Panel.from_dict(panel_dict))

from dxtbx.model import Experiment
from simtbx.nanoBragg import make_imageset
from cctbx_project.simtbx.diffBragg.phil import phil_scope
from simtbx.diffBragg import refine_launcher
E = Experiment()
E.detector = perturb_detector
E.beam = SIM.beam.nanoBragg_constructor_beam
E.crystal = C
E.imageset = make_imageset([img], E.beam, E.detector)

refls = utils.refls_from_sims([spots], E.detector, E.beam, thresh=20)

P = phil_scope.extract()
P.roi.shoebox_size = 20
P.roi.reject_edge_reflections = False
P.refiner.refine_panelRotO = [1]
P.refiner.refine_panelRotF = [1]
P.refiner.refine_panelRotS = [1]
P.refiner.refine_panelXY = [1]
P.refiner.refine_panelZ = [1]
P.refiner.ranges.panel_X = [-1e-3, 1e-3]  # meters
P.refiner.ranges.panel_Y = [-1e-3, 1e-3]
P.refiner.ranges.panel_Z = [-1.5e-3, 1e-3]
P.refiner.ranges.panel_rotO = [-10, 10]
P.refiner.ranges.panel_rotF = [-10, 10]
P.refiner.ranges.panel_rotS = [-10, 10]

P.refiner.sensitivity.panelRotOFS = [1, 1, 1]
P.refiner.sensitivity.panelXY = [1, 1]
P.refiner.sensitivity.panelZ = 1
P.refiner.max_calls = [200]
P.refiner.tradeps = 1
P.refiner.verbose = True
P.refiner.sigma_r = SIM.D.readout_noise_adu
P.refiner.plot.display = args.plot
P.refiner.plot.iteration_stride = 5
P.refiner.adu_per_photon = SIM.D.quantum_gain
P.simulator.init_scale = SIM.D.spot_scale
P.refiner.ncells_mask = "111"
P.simulator.beam.size_mm = SIM.beam.size_mm

RUC = refine_launcher.local_refiner_from_parameters(refls, E, P, miller_data=SIM.crystal.miller_array)
print("OK")
rotO, rotF, rotS = map(lambda x : 180*x / np.pi, RUC._get_panelRot_val(0))
x, y, z = RUC._get_panelXYZ_val(0)
o = abs(rotO + panel_rot_angO)
f = abs(rotF + panel_rot_angF)
s = abs(rotS + panel_rot_angS)
xx = abs(x+Xshift_mm/1000)
yy = abs(y+Yshift_mm/1000)
zz = abs(z+Zshift_mm/1000)


pixel_size = SIM.detector[0].get_pixel_size()[0]

print("RotO offset=%.8f deg." % o)
print("RotF offset=%.8f deg." % f)
print("RotS offset=%.8f deg." % s)
print("panX offset=%.8f pixels" % (xx*1e3/pixel_size))
print("panY offset=%.8f pixels" % (yy*1e3/pixel_size))
print("panZ offset=%.8f millimeters" % (zz*1e3))

assert o < 0.02
assert f < 0.02
assert s < 0.02
assert xx < 0.0025
assert yy < 0.0025
assert zz < 0.025

print("OK")
