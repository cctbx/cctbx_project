from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
parser.add_argument("--residuals", action='store_true')
parser.add_argument("--oversample", type=int, default=0)
parser.add_argument("--curvatures", action="store_true")
parser.add_argument("--nopolar", action="store_true")
args = parser.parse_args()

from dxtbx.model.crystal import Crystal
from IPython import embed
from cctbx import uctbx
from scitbx.matrix import sqr, rec, col
from dxtbx.model import Panel
from copy import deepcopy
import numpy as np
from scipy.spatial.transform import Rotation
from simtbx.diffBragg.refiners.global_refiner import GlobalRefiner
from simtbx.diffBragg.refiners.crystal_systems import MonoclinicManager

from simtbx.diffBragg.nanoBragg_crystal import nanoBragg_crystal
from simtbx.diffBragg.sim_data import SimData
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
SIM.instantiate_diffBragg(oversample=0)
SIM.D.progress_meter = False
SIM.D.verbose = 0 #1
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

# Simulate the perturbed image for comparison
panel_rot_ang = 0 #0.3
panel_rot_ang_rad = panel_rot_ang * np.pi / 180
Xshift_mm = 0.1
Yshift_mm = 0.09

det2 = deepcopy(SIM.detector)
panel = det2[0]
orig = list(panel.get_origin())
orig[0] += Xshift_mm
orig[1] += Yshift_mm
from dxtbx.model import Panel
pan_dict = panel.to_dict()
pan_dict["origin"] = orig
pan = Panel.from_dict(pan_dict)
det2[0] = pan

SIM.D.update_dxtbx_geoms(det2, SIM.beam.nanoBragg_constructor_beam, 0, panel_rot_ang_rad)
SIM.D.add_diffBragg_spots()
SIM._add_background()
SIM._add_noise()
# Perturbed image:
img_pet = SIM.D.raw_pixels.as_numpy_array()
SIM.D.raw_pixels *= 0
full_roi = SIM.D.region_of_interest

# spot_rois, abc_init , these are inputs to the refiner
# <><><><><><><><><><><><><><><><><><><><><><><><><>
spot_roi, tilt_abc = utils.process_simdata(spots, img, thresh=20) #, plot=args.plot)
n_spots = len(spot_roi)
n_kept = 30
np.random.seed(1)
idx = np.random.permutation(n_spots)#[:n_kept]
spot_roi = spot_roi[idx]
tilt_abc = tilt_abc[idx]
print ("I got %d spots!" % tilt_abc.shape[0])

#RUC = RefineDetdist(
#    spot_rois=spot_roi,
#    abc_init=tilt_abc,
#    img=img,
#    SimData_instance=SIM,
#    plot_images=args.plot,
#    plot_residuals=args.residuals)
#
#RUC.trad_conv = True
#RUC.refine_background_planes = False
#RUC.trad_conv_eps = 1e-5
#RUC.refine_detdist = True
#RUC.max_calls = 200
#RUC.run()
#
#
#
#print det2[0].get_origin()[2]
#print RUC.x[-3]
#
#assert abs(RUC.x[-3] - distance) < 1e-2
#
#print("OK!")
######3

nspot = len(spot_roi)

nanoBragg_rois = []  # special nanoBragg format
xrel, yrel, roi_imgs = [], [], []
xcom, ycom = [],[]
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
n_ncell_param = 1
nfcell_param = len(Hi_asu)
ngain_param = 1
ndetz_param = 1
REFINE_XY = True
REFINE_ROT = True
n_panRot_param = 1 * int(REFINE_ROT)
n_panXY_param = 2 * int(REFINE_XY)

n_global_unknowns = nucell_param + nfcell_param + ngain_param + ndetz_param + n_ncell_param + n_panRot_param + n_panXY_param
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

RUC.refine_background_planes =False
RUC.refine_Umatrix = False
RUC.refine_Bmatrix = False
RUC.refine_ncells =False
RUC.refine_crystal_scale = False
RUC.refine_Fcell = False
RUC.refine_detdist = False
RUC.refine_gain_fac = False
RUC.refine_panelRot = REFINE_ROT
RUC.refine_panelXY = REFINE_XY

RUC.ignore_line_search_failed_step_at_lower_bound = True

RUC.panel_group_from_id = {0: 0}  # panel group ID from panel id

init_rot = 0
if REFINE_ROT:
    init_rot = panel_rot_ang_rad
RUC.panelRot_init = {0: init_rot}
#RUC.panelRot_init = {0: panel_rot_ang_rad}  # panel group ID versus starting value
RUC.n_panel_rot_param = 1*int(RUC.refine_panelRot)
RUC.panelRot_sigma = 0.01

init_X_offset = init_Y_offset = 0
if REFINE_XY:
    init_X_offset = Xshift_mm
    init_Y_offset = Yshift_mm
RUC.panelX_init = {0: init_X_offset}
RUC.panelY_init = {0: init_Y_offset}
RUC.panelX_sigma = 0.1
RUC.panelY_sigma = 0.1
RUC.n_panel_XY_param = 2 * int(REFINE_XY)

RUC.ucell_sigmas = [1]*len(UcellMan.variables)
RUC.ucell_inits = {0:  UcellMan.variables}

RUC.max_calls = 300
RUC.trad_conv_eps = 1e-5
RUC.trad_conv = True
RUC.trial_id = 0

RUC.plot_images = args.plot
RUC.plot_stride = 1
RUC.plot_spot_stride = 1
RUC.setup_plots()

RUC.rescale_params = True
RUC.refine_rotZ = False
RUC.request_diag_once = False
RUC.S = SIM
RUC.has_pre_cached_roi_data = True
RUC.S.D.update_oversample_during_refinement = False
RUC.S.D.nopolar = False  # True

RUC.use_curvatures = False
RUC.use_curvatures_threshold = 10
RUC.calc_curvatures = args.curvatures
RUC.poisson_only = True
RUC.verbose = True
RUC.big_dump = True

RUC.run(setup_only=False)
#
#refined_distance = RUC._get_originZ_val(0)
#assert abs(refined_distance - distance) < 1e-2

#SIM.D.update_dxtbx_geoms(SIM.detector, SIM.beam.nanoBragg_constructor_beam, 0, RUC._get_panelRot_val(0))
#RUC.update
#SIM.D.update_dxtbx_geoms(SIM.detector, SIM.beam.nanoBragg_constructor_beam, 0, RUC._get_panelRot_val(0))
#SIM.D.raw_pixels *= 0
#SIM.D.region_of_interest = full_roi
#SIM.D.add_diffBragg_spots()
#SIM._add_background()
#SIM._add_noise()
## Perturbed image:
#img_opt = SIM.D.raw_pixels.as_numpy_array()
#from itertools import cycle
#imgs = cycle([img, img_opt])
#import pylab as plt
#plt.imshow(next(imgs), vmax=100)
#while 1:
#    plt.draw()
#    plt.pause(0.5)
#    plt.cla()
#    plt.imshow(next(imgs), vmax=100)
#    plt.xlim(0,200)
#    plt.ylim(200,0)

assert RUC._get_panelRot_val(0)*180 / np.pi < 0.0015
x, y = RUC._get_panelXY_val(0)

assert abs(x) < 0.0025
assert abs(y) < 0.0025

print("I AM ZIM")
print("OK!")
