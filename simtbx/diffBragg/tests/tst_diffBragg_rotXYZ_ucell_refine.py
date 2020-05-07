from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
parser.add_argument("--curvatures", action='store_true')
parser.add_argument("--rescale", action="store_true")
parser.add_argument("--poisson", action="store_true")
parser.add_argument("--readout", type=float, default=0)
parser.add_argument("--xrefinedonly", action="store_true" )
parser.add_argument("--skipfirst", action="store_true" )
parser.add_argument("--skipsecond", action="store_true" )
args = parser.parse_args()

from dxtbx.model.crystal import Crystal
from cctbx import uctbx
from IPython import embed
from scitbx.matrix import sqr, rec, col
import numpy as np
from scipy.spatial.transform import Rotation

from simtbx.diffBragg.nanoBragg_crystal import nanoBragg_crystal
from simtbx.diffBragg.sim_data import SimData
from simtbx.diffBragg import utils
from simtbx.diffBragg.refiners import RefineMissetAndUcell
from simtbx.diffBragg.refiners.crystal_systems import MonoclinicManager

ucell = (55, 65, 75, 90, 95, 90)
ucell2 = (55.1, 65.2, 74.9, 90, 94.9, 90)
#ucell2 = (55.05, 65.05, 74.05, 90, 95.05, 90)
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
nbcryst = nanoBragg_crystal()
nbcryst.dxtbx_crystal = C   # simulate ground truth
nbcryst.thick_mm = 0.1
nbcryst.Ncells_abc = 12, 12, 12

SIM = SimData()
SIM.detector = SimData.simple_detector(150, 0.1, (513, 512))
SIM.crystal = nbcryst
SIM.instantiate_diffBragg(oversample=0)
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

# Simulate the perturbed image for comparison
SIM.D.Bmatrix = C2.get_B()
SIM.D.Umatrix = C2.get_U()
SIM.D.add_diffBragg_spots()
SIM._add_background()
SIM._add_noise()
# Perturbed image:
img_pet = SIM.D.raw_pixels.as_numpy_array()
SIM.D.raw_pixels *= 0

# spot_rois, abc_init , these are inputs to the refiner
# <><><><><><><><><><><><><><><><><><><><><><><><><>
spot_roi, tilt_abc = utils.process_simdata(spots, img, thresh=20, plot=args.plot) #, edge_reflections=False)

nslow, nfast = img.shape
for i, (_, x2, _, y2) in enumerate(spot_roi):
    if x2 == nfast:
        spot_roi[i][1] = x2 - 1  # update roi_xmax
    if y2 == nslow:
        spot_roi[i][3] = y2 - 1  # update roi_ymax

UcellMan = MonoclinicManager(
    a=ucell2[0],
    b=ucell2[1],
    c=ucell2[2],
    beta=ucell2[4]*np.pi/180.)

nbcryst.dxtbx_crystal = C2
from copy import deepcopy
orig_perturbed_C2 = deepcopy(C2)
SIM.crystal = nbcryst

init_Umat_norm = np.abs(np.array(C2.get_U()) - np.array(C.get_U())).sum()
init_Bmat_norm = np.abs(np.array(C2.get_B()) - np.array(C.get_B())).sum()

if not args.skipfirst:

    RUC = RefineMissetAndUcell(
        spot_rois=spot_roi,
        abc_init=tilt_abc,
        img=img,
        SimData_instance=SIM,
        plot_images=args.plot,
        ucell_manager=UcellMan)
    RUC.trad_conv = True
    RUC.refine_background_planes = False
    RUC.refine_Amatrix = True
    RUC.trad_conv_eps = 1e-7
    RUC.max_calls = 2000
    RUC.use_curvatures = args.curvatures
    RUC.run()
    X = np.array(RUC.all_x)
    G = np.array(RUC.all_g)
    CURV = np.array(RUC.all_c)
    F = np.array(RUC.funcs)
    np.save("old", [X,G,CURV])
    np.save("oldF", F)

    ang, ax = RUC.get_correction_misset(as_axis_angle_deg=True)
    C2.rotate_around_origin(ax,ang)
    C2.set_B(RUC.get_refined_Bmatrix())

    atru, btru, ctru = C.get_real_space_vectors()
    from simtbx.diffBragg.utils import compare_with_ground_truth
    ang_off = compare_with_ground_truth(atru, btru, ctru,
                                        [C2],
                                        symbol=symbol)[0]
    print("ANG OFF= %f deg" % ang_off)


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

    # error in initial unit cell parameters
    err_init = np.linalg.norm([abs(u-u_init)/u for u, u_init in zip(ucell, ucell2)])*100

    # error in refined unit cell parameters
    err_ref = np.linalg.norm([abs(u-u_ref)/u for u, u_ref in zip(ucell, ucell_ref)])*100

    assert err_ref < 1e-1 * err_init

    # the initial perturbation matrix:
    R1 = rec(perturb_rot_axis, (3, 1)).axis_and_angle_as_r3_rotation_matrix(perturb_rot_ang, deg=True)
    # restoring purturbation applied after refinement:
    R2 = ax.axis_and_angle_as_r3_rotation_matrix(ang, deg=True)

    # the hope is that R2 cancels the effect of R1
    # hence, the product R1 and R2 should be ~ Identity
    I = np.reshape(R1*R2, (3, 3))
    assert np.all(np.round(I-np.eye(3), 3) == np.zeros((3, 3)))
    assert final_Umat_norm < 1e-1*init_Umat_norm

    print("OK")

#######################################################
if args.skipsecond:
    exit()
N_SHOTS = 1

nanoBragg_rois = []  # special nanoBragg format
xrel, yrel, roi_imgs = [], [], []
xcom, ycom = [], []
for i_roi, (x1, x2, y1, y2) in enumerate(spot_roi):
    nanoBragg_rois.append(((int(x1), int(x2)), (int(y1), int(y2))))
    yr, xr = np.indices((y2 - y1 + 1, x2 - x1 + 1))
    xrel.append(xr)
    yrel.append(yr)
    roi_imgs.append(img[y1:y2 + 1, x1:x2 + 1])
    xcom.append(.5 * (x1 + x2))
    ycom.append(.5 * (y1 + y2))

q_spot = utils.x_y_to_q(xcom, ycom, SIM.detector, SIM.beam.nanoBragg_constructor_beam)
reso = 1/np.linalg.norm(q_spot, axis=1)
all_reso = list(reso)
Ai = sqr(SIM.crystal.dxtbx_crystal.get_A()).inverse()
Ai = Ai.as_numpy_array()
HKL = np.dot(Ai, q_spot.T)
HKLi = [np.ceil(h - 0.5).astype(int) for h in HKL]
HKLi = [tuple(x) for x in np.vstack(HKLi).T]
Hi_asu = utils.map_hkl_list(HKLi, anomalous_flag=True, symbol=symbol)

#Hi_asu = [(int(h), int(k), int(l)) for h,k,l in Hi_asu]

UcellMan = MonoclinicManager(
    a=ucell2[0],
    b=ucell2[1],
    c=ucell2[2],
    beta=ucell2[4]*np.pi/180.)

nspot = len(reso)
shot_ucell_managers = {0:UcellMan}
shot_rois = {0:spot_roi}
shot_nanoBragg_rois = {0:nanoBragg_rois}
shot_roi_imgs = {0:roi_imgs}
shot_spectra = {0:SIM.beam.spectrum}
shot_crystal_GTs = {0:C}
shot_crystal_models = {0:orig_perturbed_C2}
shot_xrel = {0:xrel}
shot_yrel = {0:yrel}
shot_abc_inits = {0:tilt_abc}
shot_asu = {0:Hi_asu}  # TODO Im weird fix me
shot_hkl = {0:HKLi}  # TODO Im weird fix me
shot_panel_ids = {0: [0] * nspot}

Hi_all_ranks, Hi_asu_all_ranks = [], []
for i in range(N_SHOTS):
    Hi_all_ranks += shot_hkl[i]
    Hi_asu_all_ranks += shot_asu[i]

print("Overall completeness\n<><><><><><><><>")
from cctbx.crystal import symmetry
from cctbx import miller
uc = shot_ucell_managers[0]
from cctbx.array_family import flex as cctbx_flex
params = uc.a, uc.b, uc.c, uc.al * 180 / np.pi, uc.be * 180 / np.pi, uc.ga * 180 / np.pi
symm = symmetry(unit_cell=params, space_group_symbol=symbol)
hi_flex_unique = cctbx_flex.miller_index(list(set(Hi_asu_all_ranks)))
mset = miller.set(symm, hi_flex_unique, anomalous_flag=True)
mset.setup_binner(d_min=2, d_max=999, n_bins=10)
mset.completeness(use_binning=True).show()
print("total miller vars=%d" % (len(set(Hi_asu_all_ranks))))

# this will map the measured miller indices to their index in the LBFGS parameter array self.x
idx_from_asu = {h: i for i, h in enumerate(set(Hi_asu_all_ranks))}
# we will need the inverse map during refinement to update the miller array in diffBragg, so we cache it here
asu_from_idx = {i: h for i, h in enumerate(set(Hi_asu_all_ranks))}

# always local parameters: rotations, spot scales, tilt coeffs
nrotation_param = 3*N_SHOTS
nscale_param = 1*N_SHOTS
ntilt_param = 0  # note: tilt means tilt plane
for i_shot in range(N_SHOTS):
    ntilt_param += 3 * nspot   #nspot_per_shot[i_shot]

# unit cell parameters
nucell_param = len(shot_ucell_managers[0].variables)
n_pershot_ucell_param = nucell_param*N_SHOTS
n_global_ucell_param = 0

# mosaic domain parameter m
n_ncell_param = 1
n_pershot_m_param = 1*N_SHOTS
n_global_m_param = 0

ndetz_param = 1
n_local_unknowns = nrotation_param + nscale_param + ntilt_param + ndetz_param + n_pershot_ucell_param + n_pershot_m_param

nfcell_param = len(idx_from_asu)
ngain_param = 1

n_global_unknowns = nfcell_param + ngain_param + n_global_m_param + n_global_ucell_param
n_total_unknowns = n_local_unknowns + n_global_unknowns

from simtbx.diffBragg.refiners.global_refiner import GlobalRefiner

shot_originZ_init = {0: SIM.detector[0].get_origin()[2]}

RUC = GlobalRefiner(
    n_total_params=n_total_unknowns,
    n_local_params=n_local_unknowns,
    n_global_params=n_global_unknowns,
    local_idx_start=0,
    shot_ucell_managers=shot_ucell_managers,
    shot_rois=shot_roi_imgs,
    shot_nanoBragg_rois=shot_nanoBragg_rois,
    shot_roi_imgs=shot_roi_imgs,
    shot_spectra=shot_spectra,
    shot_crystal_GTs=shot_crystal_GTs,
    shot_crystal_models=shot_crystal_models,
    shot_xrel=shot_xrel,
    shot_yrel=shot_yrel,
    shot_abc_inits=shot_abc_inits,
    shot_asu=shot_asu,
    global_param_idx_start=n_local_unknowns,
    shot_panel_ids=shot_panel_ids,
    log_of_init_crystal_scales=None,
    all_crystal_scales=None,
    global_ncells=False,
    global_ucell=False,
    global_originZ=False,
    shot_originZ_init=shot_originZ_init,
    sgsymbol=symbol,
    omega_kahn=None)

RUC.idx_from_asu = idx_from_asu
RUC.asu_from_idx = asu_from_idx
RUC.refine_Umatrix = True
RUC.refine_Bmatrix = True
RUC.refine_ncells = False
RUC.refine_crystal_scale = False
RUC.refine_Fcell = False
RUC.refine_detdist = False
RUC.rescale_params = args.rescale
RUC.gtol = .9
RUC.trad_conv = True #False
RUC.drop_conv_max_eps = 1e-9

#RUC.ucell_sigmas = [.01, .01, .01, .01]
RUC.ucell_sigmas = 1,1,10,.1
#RUC.rotX_sigma = .001
#RUC.rotY_sigma = .001
#RUC.rotZ_sigma = .001
RUC.rotX_sigma = .1
RUC.rotY_sigma = .1
RUC.rotZ_sigma = .05
RUC.bg_coef_sigma = 1
RUC.originZ_sigma = 1
RUC.m_sigma = 1
RUC.spot_scale_sigma = 1
RUC.a_sigma = 1
RUC.b_sigma = 1
RUC.c_sigma = 1
RUC.fcell_sigma_scale = 1

RUC.max_calls = 1000
RUC.trad_conv_eps = 1e-10
RUC.trial_id = 0

RUC.plot_stride = 4
RUC.plot_spot_stride = 10
RUC.plot_residuals = False
RUC.plot_images = args.plot
RUC.setup_plots()

RUC.refine_rotZ = True
RUC.request_diag_once = False
RUC.S = SIM
if not args.curvatures:
    RUC.S.D.compute_curvatures=False
RUC.has_pre_cached_roi_data = True
RUC.S.D.update_oversample_during_refinement = False
RUC.use_curvatures = args.curvatures #False
RUC.use_curvatures_threshold = 0
RUC.calc_curvatures = True
RUC.poisson_only = False
RUC.verbose = True
RUC.big_dump = False
RUC.gt_ncells = nbcryst.Ncells_abc[0]
RUC.originZ_gt = {0:shot_originZ_init[0]}
RUC.gt_ucell = ucell[0], ucell[2]
RUC.testing_mode = False
RUC.poisson_only = args.poisson
RUC.sigma_r = args.readout #3/28.

RUC.spot_scale_init = [1]*N_SHOTS
RUC.m_init = {i: nbcryst.Ncells_abc[0] for i in range(N_SHOTS)}
RUC.ucell_inits = {i:shot_ucell_managers[i].variables for i in range(N_SHOTS)}
RUC.only_pass_refined_x_to_lbfgs = args.xrefinedonly


RUC.run(setup_only=False)
#RUC.run(setup_only=True)
if RUC.hit_break_to_use_curvatures:
    RUC.num_positive_curvatures = 0
    RUC.use_curvatures = True
    RUC.run(setup=False)

#X = np.array(RUC.all_x)
#G = np.array(RUC.all_g)
#CURV = np.array(RUC.all_c)[:,RUC.is_being_refined]
#F = np.array(RUC.funcs)
#np.save("new", [X,G,CURV])
#np.save("newF", F)
from IPython import embed
embed()

