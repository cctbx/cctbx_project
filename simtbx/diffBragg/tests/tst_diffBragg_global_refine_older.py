from __future__ import division
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
parser.add_argument("--detdist", action='store_true', help='perturb then refine the detdist')
parser.add_argument("--ncells", action='store_true', help='perturb then refine the ncells')
parser.add_argument("--spotscale", action='store_true', help='refine the crystal spot scale')
parser.add_argument("--bmatrix", action='store_true')
parser.add_argument("--umatrix", action='store_true')
parser.add_argument("--curvatures", action='store_true')
parser.add_argument("--bg", action="store_true", help='refine background planes')
parser.add_argument("--psf", action='store_true')
parser.add_argument("--gain", action='store_true')
parser.add_argument("--rescale", action="store_true")
args = parser.parse_args()

if args.detdist:
    raise NotImplementedError("Not yet refining detdist for single shots")

from dxtbx.model.crystal import Crystal
from copy import deepcopy

from dxtbx.model import Panel
from cctbx import uctbx
from scitbx.matrix import sqr, rec, col
import numpy as np
from scipy.spatial.transform import Rotation

from simtbx.diffBragg.refiners.local_refiner import LocalRefiner
from simtbx.nanoBragg import shapetype

from simtbx.nanoBragg.nanoBragg_crystal import NBcrystal
from simtbx.nanoBragg.sim_data import SimData
from simtbx.diffBragg import utils
from simtbx.diffBragg.utils import fcalc_from_pdb
from simtbx.diffBragg.refiners.crystal_systems import MonoclinicManager

ucell = (55, 65, 75, 90, 95, 90)
ucell2 = (55, 65, 75, 90, 95, 90)
if args.bmatrix:
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
perturb_rot_ang = 0
if args.umatrix:
    perturb_rot_ang = .05  # degrees

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
Ncells_gt = 12, 12, 12
nbcryst.Ncells_abc = Ncells_gt  # ground truth Ncells
nbcryst.miller_array = fcalc_from_pdb(resolution=2, wavelength=1,algorithm='fft', ucell=ucell, symbol=symbol)
print("Ground truth ncells = %f" % (nbcryst.Ncells_abc[0]))

SIM = SimData(use_default_crystal=True)
SIM.detector = SimData.simple_detector(150, 0.177, (600, 600))

# TODO get the detector model
node = SIM.detector[0]
node_d = node.to_dict()
Origin = node_d["origin"][0], node_d["origin"][1], node_d["origin"][2]
distance = Origin[2]
print ("Ground truth originZ=%f" % (SIM.detector[0].get_origin()[2]))

# TODO perturb the detector model
# copy the detector and update the origin
det2 = deepcopy(SIM.detector)
# alter the detector distance by 2 mm
detz_offset = 0.3
node_d["origin"] = Origin[0], Origin[1], Origin[2] + detz_offset
det2[0] = Panel.from_dict(node_d)

SIM.crystal = nbcryst
SIM.instantiate_diffBragg(oversample=0, auto_set_spotscale=True)
SIM.D.nopolar = False
SIM.D.default_F = 0
SIM.D.progress_meter = False
SIM.water_path_mm = 0.005
SIM.air_path_mm = 0.1
SIM.add_air = True
SIM.add_Water = True
SIM.include_noise = True
SIM.D.add_diffBragg_spots()
spots = SIM.D.raw_pixels.as_numpy_array()
SIM.D.readout_noise_adu = 3
SIM._add_background()
SIM._add_noise()

if args.psf:
    img_pre_psf = SIM.D.raw_pixels.as_numpy_array()
    v = SIM.D.verbose
    SIM.D.verbose = 8
    SIM.D.detector_psf_kernel_radius_pixels = 0
    fwhm = 76 / 56
    radius = 3
    SIM.D.apply_psf(shapetype.Fiber, fwhm, radius)
    SIM.D.verbose = v

print("Using oversample %d" % SIM.D.oversample)

# This is the ground truth image:
img = SIM.D.raw_pixels.as_numpy_array()
SIM.D.raw_pixels *= 0

if args.psf:
    y = slice(450, 480,1)
    x = slice(650, 670, 1)
    print("PSF max discrepancy: %f" % abs(img_pre_psf[y,x]- img[y,x]).max())

# Simulate the perturbed image for comparison
# perturbed detector:
if args.detdist:
    SIM.detector = det2
    SIM.D.update_dxtbx_geoms(det2, SIM.beam.nanoBragg_constructor_beam, 0)
    print("Modified originZ=%f" % (det2[0].get_origin()[2]))
# perturbed crystal:
distance = SIM.detector[0].get_origin()[2]
SIM.D.Bmatrix = C2.get_B()
SIM.D.Umatrix = C2.get_U()
nbcryst.dxtbx_crystal = C2
if args.ncells:
    Ncells_abc2 = 14, 14, 14
    nbcryst.Ncells_abc = Ncells_abc2
    SIM.D.set_value(9, Ncells_abc2[0])
    print("Modified Ncells=%f" % Ncells_abc2[0])

SIM.crystal = nbcryst
# perturbed Ncells
SIM.D.add_diffBragg_spots()
SIM._add_background()
SIM._add_noise()

# Perturbed image:
img_pet = SIM.D.raw_pixels.as_numpy_array()
SIM.D.raw_pixels *= 0

# spot_rois, abc_init , these are inputs to the refiner
# <><><><><><><><><><><><><><><><><><><><><><><><><>
spot_refls = utils.refls_from_sims(np.array([spots]), SIM.detector, SIM.beam.nanoBragg_constructor_beam, thresh=20)
out = utils.get_roi_background_and_selection_flags(spot_refls, np.array([img]), shoebox_sz=10, reject_edge_reflections=True, use_robust_estimation=True)
spot_roi, panel_ids, tilt_abc, selection_flags, background = out

UcellMan = MonoclinicManager(
    a=ucell2[0],
    b=ucell2[1],
    c=ucell2[2],
    beta=ucell2[4]*np.pi/180.)

init_Umat_norm = np.abs(np.array(C2.get_U()) - np.array(C.get_U())).sum()
init_Bmat_norm = np.abs(np.array(C2.get_B()) - np.array(C.get_B())).sum()

if args.gain:
    img = img*1.1

# TODO: the following need to be embedded in the refiner init function..
nspot = len(spot_roi)

nanoBragg_rois = []  # special nanoBragg format
xrel, yrel, roi_imgs = [], [], []
xcom, ycom = [],[]
for i_roi, (x1, x2, y1, y2) in enumerate(spot_roi):
    nanoBragg_rois.append(((x1, x2), (y1, y2)))
    yr, xr = np.indices((y2 - y1, x2 - x1))
    xrel.append(xr)
    yrel.append(yr)
    roi_imgs.append(img[y1:y2, x1:x2])

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
n_ncell_def_param = 3
n_detz_param = 1
n_eta_param = 3
nucell_param = len(UcellMan.variables)
n_ncell_param = 1
n_perspot_param = nspot
n_sausage_param = 4
n_local_unknowns = n_ncell_param + n_eta_param + nrotation_param + nscale_param + n_ncell_def_param + n_detz_param + \
                   ntilt_param + n_perspot_param + n_sausage_param

nfcell_param = len(Hi_asu)
ngain_param = 1
nspec_param = 2
npanel_param = 6  # 3 rotation and 3 offset

#n_global_unknowns = npanel_param + nspec_param + nucell_param + nfcell_param + ngain_param +  n_ncell_param
n_global_unknowns = npanel_param + nspec_param + nfcell_param + ngain_param + nucell_param
n_total_unknowns = n_local_unknowns + n_global_unknowns

RUC = LocalRefiner(
    n_total_params=n_total_unknowns,
    n_local_params=n_local_unknowns,
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
    perturb_fcell=False,
    global_ncells=False,
    global_ucell=True,
    global_detector_distance=False,
    shot_detector_distance_init={0: 0},
    sgsymbol=symbol)

#TODO make this part of class init:
idx_from_asu = {h: i for i, h in enumerate(set(Hi_asu))}
asu_from_idx = {i: h for i, h in enumerate(set(Hi_asu))}

RUC.idx_from_asu = idx_from_asu
RUC.asu_from_idx = asu_from_idx

RUC.refine_background_planes = args.bg
RUC.refine_Umatrix = args.umatrix
RUC.refine_Bmatrix = args.bmatrix
RUC.refine_ncells = args.ncells
RUC.refine_crystal_scale = args.spotscale
RUC.refine_Fcell = False
RUC.selection_flags = {0: selection_flags}
RUC.refine_detdist = args.detdist
RUC.refine_gain_fac = args.gain
RUC.ucell_sigmas = [.1, .1, .3, .005]
RUC.ucell_inits = {0: UcellMan.variables}
RUC.m_init = {0:SIM.crystal.Ncells_abc[0]}
RUC.spot_scale_init = {0: 1}

RUC.max_calls = 300
RUC.trad_conv_eps = 1e-2
RUC.trad_conv = True
RUC.trial_id = 0

RUC.plot_stride = 10
RUC.plot_residuals = False
RUC.plot_images = args.plot
RUC.setup_plots()

RUC.refine_rotZ = True
RUC.request_diag_once = False
RUC.S = SIM
RUC.has_pre_cached_roi_data = True
RUC.S.D.update_oversample_during_refinement = False
RUC.use_curvatures = False
RUC.use_curvatures_threshold = 6
RUC.calc_curvatures = args.curvatures
RUC.poisson_only = False #True
RUC.verbose = True
RUC.big_dump = True
RUC.rescale_params = args.rescale

RUC.gt_ucell = ucell[0], ucell[1], ucell[2], ucell[4]
RUC.gt_ncells = Ncells_gt[0]
RUC.testing_mode = True

RUC.bg_offset_only = True
RUC.bg_offset_positive = True

RUC.run(setup_only=False)
if RUC.hit_break_to_use_curvatures:
    RUC.num_positive_curvatures = 0
    RUC.use_curvatures = True
    RUC.run(setup=False)

ang, ax = RUC.get_correction_misset(as_axis_angle_deg=True, i_shot=0)
if ang > 0:
    C2.rotate_around_origin(ax, ang)
C2.set_B(RUC.get_refined_Bmatrix(i_shot=0))

final_Umat_norm = np.abs(np.array(C2.get_U()) - np.array(C.get_U())).sum()
final_Bmat_norm = np.abs(np.array(C2.get_B()) - np.array(C.get_B())).sum()

# refined unit cell parameters
ucell_ref = C2.get_unit_cell().parameters()

print("Results!")
print("Before refinement: Umatrix distance=%2.7g, Bmatrix distance=%2.7g" % (init_Umat_norm, init_Bmat_norm))
print("After refinement: Umatrix distance=%2.7g, Bmatrix distance=%2.7g" % (final_Umat_norm, final_Bmat_norm))
print("")
print("ground truth unit cell: %2.7g,%2.7g,%2.7g,%2.7g,%2.7g,%2.7g" % ucell)
print("unit cell passed to refinement: %2.7g,%2.7g,%2.7g,%2.7g,%2.7g,%2.7g" % ucell2)
print("refined unit cell: %2.7g,%2.7g,%2.7g,%2.7g,%2.7g,%2.7g" % ucell_ref)
print("")
print("Perturbation axis =%+2.7g,%+2.7g,%+2.7g and angle=%+2.7g deg"
      % (perturb_rot_axis[0], perturb_rot_axis[1], perturb_rot_axis[2], perturb_rot_ang))
print("Misset applied during refinement: axis=%+2.7g,%+2.7g,%+2.7g and angle=%+2.7g deg"
      % (ax[0], ax[1], ax[2], ang))

print("OK")
