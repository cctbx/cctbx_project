from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
args = parser.parse_args()

from dxtbx.model.crystal import Crystal
from cctbx import uctbx
from scitbx.matrix import sqr, rec, col
from scipy.spatial.transform import Rotation
import numpy as np

from simtbx.diffBragg.nanoBragg_crystal import nanoBragg_crystal
from simtbx.diffBragg.sim_data import SimData
from simtbx.diffBragg import utils
from simtbx.diffBragg.refiners import RefineNcells


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
nbcryst = nanoBragg_crystal()
nbcryst.dxtbx_crystal = C   # simulate ground truth
nbcryst.thick_mm = 0.1
nbcryst.Ncells_abc = 19, 19, 19

print "Ground truth ncells abc=%f" % (nbcryst.Ncells_abc[0])

# generate the ground truth image
SIM = SimData()
SIM.detector = SimData.simple_detector(200, 0.1, (1024, 1024))
SIM.crystal = nbcryst
SIM.instantiate_diffBragg(oversample=0)
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
SIM.crystal.Ncells_abc = Ncells_abc_2
SIM.D.set_value(9, Ncells_abc_2[0])
SIM.D.add_diffBragg_spots()
SIM._add_background()
SIM._add_noise()
# Perturbed image:
img_pet = SIM.D.raw_pixels.as_numpy_array()
SIM.D.raw_pixels *= 0

# spot_rois, abc_init , these are inputs to the refiner
# <><><><><><><><><><><><><><><><><><><><><><><><><>
spot_roi, tilt_abc = utils.process_simdata(spots, img, thresh=20, plot=args.plot)
n_spots = len(spot_roi)
n_kept = 20
np.random.seed(1)
idx = np.random.permutation(n_spots)[:n_kept]
spot_roi = spot_roi[idx]
tilt_abc = tilt_abc[idx]

RUC = RefineNcells(
    spot_rois=spot_roi,
    abc_init=tilt_abc,
    img=img,
    SimData_instance=SIM,
    plot_images=args.plot,
    plot_residuals=True)

RUC.refine_Amatrix = False
RUC.refine_ncells = True
RUC.refine_gain_fac = False
RUC.refine_crystal_scale = False
RUC.trad_conv = True
RUC.trad_conv_eps = 1e-5
RUC.max_calls = 100
RUC.run()

assert round(np.exp(RUC.x[-3])) == 19
print("OK.")
