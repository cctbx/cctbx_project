from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
parser.add_argument("--curvatures", action='store_true')
args = parser.parse_args()

from dxtbx.model.crystal import Crystal
from cctbx import uctbx
from scitbx.matrix import sqr, rec, col
from scipy.spatial.transform import Rotation
import numpy as np

from simtbx.diffBragg.nanoBragg_crystal import nanoBragg_crystal
from simtbx.diffBragg.sim_data import SimData
from simtbx.diffBragg import utils
from simtbx.diffBragg.refiners.crystal_systems import MonoclinicManager
from simtbx.diffBragg.refiners.global_refiner import GlobalRefiner


ucell = (85.2, 96, 124, 90, 105, 90)
symbol = "P121"

UcellMan = MonoclinicManager(
    a=ucell[0],
    b=ucell[1],
    c=ucell[2],
    beta=ucell[4]*np.pi/180.)

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
nbcryst.isotropic_ncells = False
ncells_gt = 19, 23, 15
nbcryst.Ncells_abc = ncells_gt

print("Ground truth ncells abc=%f" % (nbcryst.Ncells_abc[0]))

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
print("I got %d spots!" % tilt_abc.shape[0])

nspot = len(spot_roi)

nanoBragg_rois = []  # special nanoBragg format
xrel, yrel, roi_imgs = [], [], []
xcom, ycom = [], []
for i_roi, (x1, x2, y1, y2) in enumerate(spot_roi):
    nanoBragg_rois.append(((x1, x2), (y1, y2)))
    yr, xr = np.indices((y2 - y1 + 1, x2 - x1 + 1))
    xrel.append(xr)
    yrel.append(yr)
    roi_imgs.append(img[y1:y2 + 1, x1:x2 + 1])

    xcom.append(.5*(x1 + x2))
    ycom.append(.5*(x1 + x2))

q_spot = utils.x_y_to_q(xcom, ycom, SIM.detector, SIM.beam.nanoBragg_constructor_beam)
Ai = sqr(SIM.crystal.dxtbx_crystal.get_A()).inverse()
Ai = Ai.as_numpy_array()
HKL = np.dot(Ai, q_spot.T)
HKLi = [np.ceil(h - 0.5).astype(int) for h in HKL]
HKLi = [tuple(x) for x in np.vstack(HKLi).T]
Hi_asu = utils.map_hkl_list(HKLi, anomalous_flag=True, symbol=symbol)

nrotation_param = 3
nscale_param = 1
ntilt_param = 3*nspot
n_local_unknowns = nrotation_param + nscale_param + ntilt_param

UcellMan = MonoclinicManager(a=ucell[0], b=ucell[1], c=ucell[2], beta=ucell[4]*np.pi/180.)
nucell_param = len(UcellMan.variables)
n_ncell_param = 3
nfcell_param = len(Hi_asu)
ngain_param = 1
ndetz_param = 1

n_global_unknowns = nucell_param + nfcell_param + ngain_param + ndetz_param + n_ncell_param
n_total_unknowns = n_local_unknowns + n_global_unknowns

SIM.D.oversample_omega = False
starting_originZ = SIM.detector[0].get_origin()[2]
RUC = GlobalRefiner(
    n_total_params=n_total_unknowns,
    n_local_params=n_local_unknowns,
    n_global_params=n_global_unknowns,
    local_idx_start=0,
    shot_ucell_managers={0: UcellMan},
    shot_rois={0: spot_roi},
    shot_nanoBragg_rois={0: nanoBragg_rois},
    shot_roi_imgs={0: roi_imgs},
    shot_spectra={0: SIM.beam.spectrum},
    shot_crystal_GTs={0: C},
    shot_crystal_models={0: SIM.crystal.dxtbx_crystal},
    shot_xrel={0: xrel},
    shot_yrel={0: yrel},
    shot_abc_inits={0: tilt_abc},
    shot_asu={0: Hi_asu},
    global_param_idx_start=n_local_unknowns,
    shot_panel_ids={0: [0]*nspot},
    log_of_init_crystal_scales=None,
    all_crystal_scales=None,
    shot_originZ_init={0: starting_originZ},
    perturb_fcell=False,
    global_ncells=True,
    global_ucell=True,
    sgsymbol=symbol)

#TODO make this part of class init:
idx_from_asu = {h: i for i, h in enumerate(set(Hi_asu))}
asu_from_idx = {i: h for i, h in enumerate(set(Hi_asu))}

RUC.idx_from_asu = idx_from_asu
RUC.asu_from_idx = asu_from_idx

RUC.ucell_inits = {0:  UcellMan.variables}
RUC.ucell_sigmas = [1 for _ in UcellMan.variables]
RUC.refine_background_planes = False
RUC.refine_Umatrix = False
RUC.refine_Bmatrix = False
RUC.refine_ncells = True
RUC.refine_crystal_scale = False
RUC.refine_Fcell = False
RUC.refine_detdist = False
RUC.refine_gain_fac = False
RUC.n_ncells_param = 3  # TODO : make this implicit somehow?
RUC.m_init = {0: Ncells_abc_2}  # either a float or a 3-tuple of floats

RUC.max_calls = 300
RUC.trad_conv_eps = 1e-5
RUC.trad_conv = True
RUC.trial_id = 0

RUC.plot_images = False
RUC.setup_plots()

RUC.rescale_params = True
RUC.refine_rotZ = False
RUC.request_diag_once = False
RUC.S = SIM
RUC.has_pre_cached_roi_data = True
RUC.S.D.update_oversample_during_refinement = False
RUC.use_curvatures = False
RUC.use_curvatures_threshold = 10
RUC.calc_curvatures = args.curvatures
RUC.poisson_only = True
RUC.verbose = True
RUC.big_dump = True

RUC.run(setup_only=False)
if RUC.hit_break_to_use_curvatures:
    RUC.num_positive_curvatures = 0
    RUC.use_curvatures = True
    RUC.run(setup=False)
Ncells_opt = RUC._get_m_val(0)
print("Ncells optimized: ", Ncells_opt)
print("Ncells GT: ", ncells_gt)

print("I AM ZIM")
print("OK!")

