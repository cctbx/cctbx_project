from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
args = parser.parse_args()

from dxtbx.model.crystal import Crystal
from IPython import embed
from cctbx import uctbx
from scitbx.matrix import sqr, rec, col
from dxtbx.model import Panel
from copy import deepcopy
import numpy as np
from scipy.spatial.transform import Rotation
import pylab as plt

from simtbx.diffBragg.nanoBragg_crystal import nanoBragg_crystal
from simtbx.diffBragg.sim_data import SimData
from simtbx.diffBragg import utils
from simtbx.diffBragg.refiners import RefineDetdist


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
nbcryst.Ncells_abc = 12, 12, 12

SIM = SimData()
SIM.detector = SimData.simple_detector(160, 0.1, (1024, 1024))

# grab the detector node (ground truth)
node = SIM.detector[0]
node_d = node.to_dict()
Origin = node_d["origin"][0], node_d["origin"][1], node_d["origin"][2]
distance = Origin[2]
print "Ground truth originZ=%f" % (SIM.detector[0].get_origin()[2])


# copy the detector and update the origin
det2 = deepcopy(SIM.detector)
# alter the detector distance by 2 mm
node_d["origin"] = Origin[0], Origin[1], Origin[2]+3
det2[0] = Panel.from_dict(node_d)
print ("Modified originZ=%f" % (det2[0].get_origin()[2]))

SIM.crystal = nbcryst
SIM.instantiate_diffBragg(oversample=0)
SIM.D.progress_meter = False
SIM.D.verbose = 0 #1
SIM.D.nopolar = True
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

# Simulate the perturbed image for comparison
SIM.detector = det2
SIM.D.update_dxtbx_geoms(det2, SIM.beam.nanoBragg_constructor_beam, 0)
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
import numpy as np
n_kept = 20
np.random.seed(1)
idx = np.random.permutation(n_spots)[:n_kept]
spot_roi = spot_roi[idx]
tilt_abc = tilt_abc[idx]
print ("I got %d spots!" % tilt_abc.shape[0])

RUC = RefineDetdist(
    spot_rois=spot_roi,
    abc_init=tilt_abc,
    img=img,
    SimData_instance=SIM,
    plot_images=args.plot,
    plot_residuals=True)

RUC.trad_conv = True
RUC.refine_background_planes = False
RUC.trad_conv_eps = 1e-5
RUC.max_calls = 200
RUC._setup()
RUC._cache_roi_arrays()
#RUC.run()

from scitbx.array_family import flex

def func(x, RUC):
    RUC.x = flex.double(x)
    f, g = RUC.compute_functional_and_gradients()
    return f


def fprime(x, RUC):
    RUC.x = flex.double(x)
    f, g = RUC.compute_functional_and_gradients()
    return g.as_numpy_array()


from scipy.optimize import fmin_l_bfgs_b

bounds = [(-np.inf, np.inf)]*RUC.n
bounds[-3] = -170, -140


print("GO!")
out = fmin_l_bfgs_b(func=func, x0=np.array(RUC.x), fprime=fprime, args=[RUC], bounds=bounds)
from IPython import embed
embed()

