from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
parser.add_argument("--detdist", action='store_true', help='perturb then refine the detdist')
parser.add_argument("--ncells", action='store_true', help='perturb then refine the ncells')
parser.add_argument("--bg", action='store_true', help='refine bg planes... ')
parser.add_argument("--spotscale", action='store_true')
parser.add_argument("--poisson", action='store_true')
parser.add_argument("--bmatrix", action='store_true')
parser.add_argument("--umatrix", action='store_true')
parser.add_argument("--fcell", action='store_true')
parser.add_argument("--nshots", default=1, type=int)
parser.add_argument("--curvatures", action='store_true')
parser.add_argument("--psf", action='store_true')
parser.add_argument("--gain", action='store_true')
parser.add_argument("--iterfreeze", action="store_true")
parser.add_argument("--rescale", action="store_true")
parser.add_argument("--onlyindexed", action="store_true")
parser.add_argument("--testbg", action="store_true")
parser.add_argument("--testfcell", action="store_true")
parser.add_argument("--testUmatrix", action="store_true")
parser.add_argument("--bgoffsetonly", action="store_true")
parser.add_argument("--bgoffsetpositive", action="store_true")
parser.add_argument("--shufflebg", action="store_true")
parser.add_argument("--tiltfit", action="store_true")
parser.add_argument("--predictwithtruth", action="store_true")
parser.add_argument("--pershotucell", action="store_true")
parser.add_argument("--pershotncells", action="store_true")
parser.add_argument("--maxcalls", type=int, default=None)
parser.add_argument("--displayedimage", type=int, default=0)
parser.add_argument("--perturbfcell", type=float, default=None, help="perturbation factor, 0.1 is small, 1 is large")
parser.add_argument("--fractionperturbed", type=float, default=0.1, help="Fraction of Fhkl to perturn")
parser.add_argument("--fcellsigmascale", type=float, default=None)
args = parser.parse_args()


from dxtbx.model.crystal import Crystal
from copy import deepcopy

from dxtbx.model import Panel
from cctbx import uctbx
from scitbx.matrix import sqr, rec, col
import numpy as np
from scipy.spatial.transform import Rotation
from scitbx.matrix import sqr

from simtbx.nanoBragg import shapetype
from simtbx.diffBragg.nanoBragg_crystal import nanoBragg_crystal
from simtbx.diffBragg.sim_data import SimData
from simtbx.diffBragg import utils
from simtbx.diffBragg.refiners.global_refiner import GlobalRefiner
from IPython import embed
from simtbx.diffBragg.refiners.crystal_systems import MonoclinicManager, TetragonalManager

# containers for GlobalRefine
shot_ucell_managers={}
shot_rois={}
shot_nanoBragg_rois={}
shot_roi_imgs={}
shot_spectra={}
shot_crystal_GTs={}
shot_crystal_models={}
shot_xrel={}
shot_yrel={}
shot_abc_inits={}
shot_asu={}
shot_hkl={}
shot_panel_ids={}
nspot_per_shot = {}
shot_originZ_init = {}
# GLOBAL PARAMETERS

all_c_before = []

ucell = (79, 79, 38, 90,90,90)
symbol = "P43212"

from simtbx.diffBragg.utils import fcalc_from_pdb
miller_array_GT = fcalc_from_pdb(resolution=2, wavelength=1, algorithm='fft', ucell=ucell, symbol=symbol)
Ncells_gt = 12, 12, 12

N_SHOTS = args.nshots

np.random.seed(3142019)
detdists_gt = np.random.normal(150,0.1, N_SHOTS)
offsets = np.random.uniform(1, 3, N_SHOTS) * np.random.choice([1,-1], N_SHOTS)
originZ_gt = {}
all_dets = []
all_reso = []

for i_shot in range(N_SHOTS):

    # generate a random raotation
    rotation = Rotation.random(num=1, random_state=100 + i_shot)[0]
    Q = rec(rotation.as_quat(), n=(4, 1))
    rot_ang, rot_axis = Q.unit_quaternion_as_axis_and_angle()

    # generate a small perturbation rotation
    perturb_rot_axis = np.random.random(3)
    perturb_rot_axis /= np.linalg.norm(perturb_rot_axis)
    perturb_rot_ang = 0
    #if args.umatrix:
    #    perturb_rot_ang = np.random.choice([0.02, 0.03, 0.04, .05])  # degrees

    # make the ground truth crystal:
    a_real, b_real, c_real = sqr(uctbx.unit_cell(ucell).orthogonalization_matrix()).transpose().as_list_of_lists()
    x = col((-1, 0, 0))
    y = col((0, -1, 0))
    z = col((0, 0, -1))
    rx, ry, rz = np.random.uniform(-180, 180, 3)
    RX = x.axis_and_angle_as_r3_rotation_matrix(rx, deg=True)
    RY = y.axis_and_angle_as_r3_rotation_matrix(ry, deg=True)
    RZ = z.axis_and_angle_as_r3_rotation_matrix(rz, deg=True)
    M = RX * RY * RZ
    a_real = M * col(a_real)
    b_real = M * col(b_real)
    c_real = M * col(c_real)
    C = Crystal(a_real, b_real, c_real, symbol)
    C.rotate_around_origin(rot_axis, rot_ang)

    # Setup the simulation and create a realistic image
    # with background and noise
    # <><><><><><><><><><><><><><><><><><><><><><><><><>
    nbcryst = nanoBragg_crystal()
    nbcryst.dxtbx_crystal = C   # simulate ground truth
    nbcryst.thick_mm = 0.1
    nbcryst.Ncells_abc = Ncells_gt  # ground truth Ncells

    nbcryst.miller_array = miller_array_GT
    print("Ground truth ncells = %f" % (nbcryst.Ncells_abc[0]))

    # ground truth detector
    DET_gt = SimData.simple_detector(detdists_gt[i_shot], 0.177, (600, 600))
    originZ_gt[i_shot] = DET_gt[0].get_origin()[2]

    # initialize the simulator
    SIM = SimData()
    SIM.detector = DET_gt
    all_dets.append(DET_gt)

    node = SIM.detector[0]
    node_d = node.to_dict()
    Origin = node_d["origin"][0], node_d["origin"][1], node_d["origin"][2]
    distance = Origin[2]
    print("Ground truth originZ=%f" % (SIM.detector[0].get_origin()[2]))
    shot_originZ_init[i_shot] = distance

    SIM.crystal = nbcryst
    SIM.instantiate_diffBragg(oversample=0)
    SIM.D.nopolar = False
    SIM.D.default_F = 0
    SIM.D.progress_meter = False
    if args.curvatures:
        SIM.D.compute_curvatures = True
    SIM.water_path_mm = 0.15
    SIM.air_path_mm = 0.1
    SIM.D.F000 = 0
    SIM.D.add_diffBragg_spots()
    SPOTS = SIM.D.raw_pixels.as_numpy_array()
    SIM.D.readout_noise_adu = 1
    SIM._add_background()
    BACKGROUND_IMAGE = SIM.D.raw_pixels.as_numpy_array() - SPOTS

    # This is the ground truth image:
    img = SIM.D.raw_pixels.as_numpy_array()
    # get the polarization and kahn facrors
    SIM.D.raw_pixels *= 0
    SIM.D.only_save_omega_kahn = True
    SIM.D.add_diffBragg_spots()
    omega_kahn = SIM.D.raw_pixels.as_numpy_array()
    SIM.D.raw_pixels *= 0
    SIM.D.only_save_omega_kahn = False

    SIM.D.Bmatrix = C.get_B()
    SIM.D.Umatrix = C.get_U()
    nbcryst.dxtbx_crystal = C

    nbcryst.Ncells_abc = Ncells_gt
    SIM.D.set_value(9, Ncells_gt[0])

    SIM.crystal = nbcryst
    SIM.D.raw_pixels *= 0
    SIM.D.add_diffBragg_spots()
    SPOTS2 = SIM.D.raw_pixels.as_numpy_array()
    SIM.D.raw_pixels *= 0

    if args.tiltfit:
        from tilt_fit.tilt_fit import tilt_fit
        from cxid9114.prediction import prediction_utils

        expLst = SIM.D.as_explist()
        exper = expLst[0]
        if args.predictwithtruth:
            refls_predict = prediction_utils.refls_from_sims([SPOTS], exper.detector, exper.beam, thresh=20)
        else:
            refls_predict = prediction_utils.refls_from_sims([SPOTS2], exper.detector, exper.beam, thresh=20)
        results = tilt_fit(
            imgs=[img/omega_kahn], is_bg_pix=[SPOTS < 20],
            delta_q=0.095, photon_gain=1, sigma_rdout=1, zinger_zscore=5,
            exper=exper, predicted_refls=refls_predict, sb_pad=2)

        refls_predict, tilt_abc, error_in_tilt, I_Leslie99, varI_Leslie99 = results
        shoeboxes = refls_predict['shoebox']
        spot_roi = np.vstack([list(sb.bbox)[:4] for sb in shoeboxes])

        Hi = np.vstack(refls_predict['miller_index'])
        did_index = np.array(refls_predict['id']) != -1  # refls that didnt index should be labeled with -1
        boundary_spot = np.array(refls_predict['boundary'])
        resolution = np.array(refls_predict["resolution"])  # reso of the spots
        if args.onlyindexed:
            spot_roi = spot_roi[did_index]
            tilt_abc = tilt_abc[did_index]
            Hi = Hi[did_index]
            resolution = resolution[did_index]

    else:
        if args.predictwithtruth:
            spot_roi, tilt_abc = utils.process_simdata(SPOTS, img/omega_kahn, thresh=20, plot=args.plot, shoebox_sz=20)
        else:
            spot_roi, tilt_abc = utils.process_simdata(SPOTS2, img/omega_kahn, thresh=20, plot=args.plot, shoebox_sz=20)
        #spot_roi, tilt_abc = utils.process_simdata(SPOTS, img, thresh=20, plot=args.plot, shoebox_sz=20)

    UcellMan = TetragonalManager(
        a=ucell[0],
        c=ucell[2])

    # TODO: the following need to be added to the refiner init function..
    nspot = len(spot_roi)
    nspot_per_shot[i_shot] = nspot

    nanoBragg_rois = []  # special nanoBragg format
    xrel, yrel, roi_imgs = [], [], []
    xcom, ycom = [],[]
    for i_roi, (x1, x2, y1, y2) in enumerate(spot_roi):
        nanoBragg_rois.append(((int(x1), int(x2)), (int(y1), int(y2))))
        yr, xr = np.indices((y2 - y1 + 1, x2 - x1 + 1))
        xrel.append(xr)
        yrel.append(yr)
        roi_imgs.append(img[y1:y2 + 1, x1:x2 + 1])
        xcom.append(.5*(x1 + x2))
        ycom.append(.5*(y1 + y2))

    if args.tiltfit:
        HKLi = [tuple(hi) for hi in Hi]
    else:
        q_spot = utils.x_y_to_q(xcom, ycom, SIM.detector, SIM.beam.nanoBragg_constructor_beam)
        reso = 1/np.linalg.norm(q_spot, axis=1)
        all_reso += list(reso)
        Ai = sqr(SIM.crystal.dxtbx_crystal.get_A()).inverse()
        Ai = Ai.as_numpy_array()
        HKL = np.dot(Ai, q_spot.T)
        HKLi = [np.ceil(h - 0.5).astype(int) for h in HKL]
        HKLi = [tuple(x) for x in np.vstack(HKLi).T]

    Hi_asu = utils.map_hkl_list(HKLi, anomalous_flag=True, symbol=symbol)

    shot_ucell_managers[i_shot]= UcellMan
    shot_rois[i_shot]= spot_roi
    shot_nanoBragg_rois[i_shot]= nanoBragg_rois
    shot_roi_imgs[i_shot]= roi_imgs
    shot_spectra[i_shot]= SIM.beam.spectrum
    shot_crystal_GTs[i_shot]= C
    shot_crystal_models[i_shot]= SIM.crystal.dxtbx_crystal
    shot_xrel[i_shot]= xrel
    shot_yrel[i_shot]= yrel
    shot_abc_inits[i_shot]= tilt_abc
    shot_asu[i_shot]= Hi_asu  # TODO Im weird fix me
    shot_hkl[i_shot]= HKLi  # TODO Im weird fix me
    shot_panel_ids[i_shot]= [0]*nspot

    if i_shot < N_SHOTS-1:
        SIM.D.free_all()

#if args.detdist:
#    SIM.D.oversample_omega = False  # necessary to refine detector distance

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
    ntilt_param += 3 * nspot_per_shot[i_shot]

# unit cell parameters
nucell_param = len(shot_ucell_managers[0].variables)
n_pershot_ucell_param = 0
n_global_ucell_param = nucell_param
if args.pershotucell:
    n_pershot_ucell_param += nucell_param*N_SHOTS
    n_global_ucell_param = 0

# mosaic domain parameter m
n_ncell_param = 1
n_pershot_m_param = 0
n_global_m_param = n_ncell_param
if args.pershotncells:
    n_pershot_m_param = 1*N_SHOTS
    n_global_m_param = 0

ndetz_param = len(detdists_gt)
n_local_unknowns = nrotation_param + nscale_param + ntilt_param + ndetz_param + n_pershot_ucell_param + n_pershot_m_param

nfcell_param = len(idx_from_asu)
ngain_param = 1

n_global_unknowns = nfcell_param + ngain_param + n_global_m_param + n_global_ucell_param
n_total_unknowns = n_local_unknowns + n_global_unknowns

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
    global_ncells=not args.pershotncells,
    global_ucell=not args.pershotucell,
    global_originZ=False,
    shot_originZ_init=shot_originZ_init,
    sgsymbol=symbol,
    omega_kahn=[omega_kahn])

RUC.idx_from_asu = idx_from_asu
RUC.asu_from_idx = asu_from_idx
RUC.refine_background_planes = args.bg
RUC.refine_Umatrix = args.umatrix
RUC.refine_Bmatrix = args.bmatrix
RUC.refine_ncells = args.ncells
RUC.refine_crystal_scale = args.spotscale
RUC.refine_Fcell = args.fcell
RUC.refine_detdist = args.detdist
RUC.refine_gain_fac = args.gain
RUC.rescale_params = args.rescale
RUC.max_calls = 1000
if args.maxcalls is not None:
    RUC.max_calls = args.maxcalls
RUC.trad_conv_eps = 1e-7
RUC.trad_conv = True
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
    RUC.S.D.compute_curvatures = False
RUC.has_pre_cached_roi_data = True
RUC.S.D.update_oversample_during_refinement = False
RUC.use_curvatures = False
RUC.use_curvatures_threshold = 10
RUC.bg_offset_positive = args.bgoffsetpositive
RUC.bg_offset_only = args.bgoffsetonly
RUC.calc_curvatures = args.curvatures
RUC.poisson_only = args.poisson
RUC.verbose = True
RUC.big_dump = False
RUC.gt_ncells = Ncells_gt[0]
RUC.originZ_gt = originZ_gt
RUC.gt_ucell = ucell[0], ucell[2]
RUC.spot_scale_init = [1]*N_SHOTS
RUC.testing_mode = False

RUC.m_init = {0:Ncells_gt[0]}  # np.log(Ncells_abc2[0]-3)
RUC.ucell_inits = {0: [ucell[0], ucell[2]]}

#RUC.S.D.update_oversample_during_refinement = False  # todo: decide
Fobs = RUC.S.crystal.miller_array_high_symmetry
RUC.Fref = miller_array_GT
#dmax, dmin = Fobs.d_max_min()
dmax, dmin = max(all_reso), min(all_reso)
RUC.binner_dmax = dmax + 1e-6
RUC.binner_dmin = dmin - 1e-6
RUC.binner_nbin = 10
RUC.scale_r1 = True
RUC.merge_stat_frequency = 1 #2
RUC.print_resolution_bins = False
if args.fcellsigmascale is not None:
    RUC.fcell_sigma_scale = args.fcellsigmascale

RUC.run(setup_only=True)

###########################################
###########################################
###########################################

#def compute_functional_and_gradients(self):

from cctbx.array_family import flex
# reset gradient and functional
RUC.target_functional = 0
RUC.grad = flex.double(RUC.n)
if RUC.calc_curvatures:
    RUC.curv = flex.double(RUC.n)

# current work has these all at 1
RUC.gain_fac = RUC.x[RUC.gain_xpos]
RUC.G2 = RUC.gain_fac ** 2

RUC._update_Fcell()  # update the structure factor with the new x

RUC.i_shot = 0

RUC.scale_fac = RUC._get_spot_scale(RUC._i_shot)

# TODO: Omatrix update? All crystal models here should have the same to_primitive operation, ideally
RUC._update_beams()
RUC._update_umatrix()
RUC._update_ucell()
RUC._update_ncells()
RUC._update_rotXYZ()
n_spots = len(RUC.NANOBRAGG_ROIS[RUC._i_shot])
i_spot = 2

RUC._panel_id = RUC.PANEL_IDS[RUC._i_shot][i_spot]

RUC.Imeas = RUC.ROI_IMGS[RUC._i_shot][i_spot]
RUC._update_dxtbx_detector()
RUC._run_diffBragg_current(i_spot)
RUC._set_background_plane(i_spot)
RUC._extract_pixel_data()
RUC._evaluate_averageI()
img0 = RUC.model_Lambda.copy()

if RUC.poisson_only:
    RUC._evaluate_log_averageI()
else:
    RUC._evaluate_log_averageI_plus_sigma_readout()

RUC._derivative_convenience_factors()
L = RUC._target_accumulate() #RUC._gaussian_target()
d, d2 = RUC._spot_scale_derivatives(return_derivatives=True)
dL_dtheta = RUC._grad_accumulate(d)
d2L_dtheta2 = RUC._curv_accumulate(d, d2)

all_error = []
all_Lerror = []
all_error2 = []
all_Lerror2 = []
shifts = [1e-6 * (2**i) for i in range(10)]
all_delta_h = []
all_delta_h2 = []
for delta_s in shifts:
    # update spot scale
    spot_scale_init = RUC.spot_scale_init[i_shot]
    new_spot_scale = spot_scale_init + delta_s
    RUC._set_spot_scale(new_spot_scale, i_shot)

    # NOTE if not args.rescale
    delta_h = np.log(new_spot_scale) - np.log(spot_scale_init)
    all_delta_h.append(delta_h)

    # redo the target function
    RUC.scale_fac = RUC._get_spot_scale(RUC._i_shot)

    # TODO: Omatrix update? All crystal models here should have the same to_primitive operation, ideally
    RUC._update_beams()
    RUC._update_umatrix()
    RUC._update_ucell()
    RUC._update_ncells()
    RUC._update_rotXYZ()
    n_spots = len(RUC.NANOBRAGG_ROIS[RUC._i_shot])
    i_spot = 2
    RUC._panel_id = RUC.PANEL_IDS[RUC._i_shot][i_spot]
    RUC.Imeas = RUC.ROI_IMGS[RUC._i_shot][i_spot]
    RUC._update_dxtbx_detector()
    RUC._run_diffBragg_current(i_spot)
    RUC._set_background_plane(i_spot)
    RUC._extract_pixel_data()
    RUC._evaluate_averageI()

    img_plus = RUC.model_Lambda.copy()
    img_fdiff = (img_plus -img0)/delta_h
    error = np.abs(d - img_fdiff).mean()

    if RUC.poisson_only:
        RUC._evaluate_log_averageI()
    else:
        RUC._evaluate_log_averageI_plus_sigma_readout()

    RUC._derivative_convenience_factors()
    Lplus = RUC._target_accumulate() #RUC._gaussian_target()

    finite_difference = (Lplus - L) / delta_h
    Lerror = abs(finite_difference - dL_dtheta)

    if args.curvatures:
        new_spot_scale = np.exp(RUC.x[RUC.spot_scale_xpos[0]] - 2*delta_h)
        RUC._set_spot_scale(new_spot_scale, i_shot)
        # NOTE if not args.rescale

        # redo the target function
        RUC.scale_fac = RUC._get_spot_scale(RUC._i_shot)
        RUC._update_beams()
        RUC._update_umatrix()
        RUC._update_ucell()
        RUC._update_ncells()
        RUC._update_rotXYZ()
        n_spots = len(RUC.NANOBRAGG_ROIS[RUC._i_shot])
        i_spot = 2
        RUC._panel_id = RUC.PANEL_IDS[RUC._i_shot][i_spot]
        RUC.Imeas = RUC.ROI_IMGS[RUC._i_shot][i_spot]
        RUC._update_dxtbx_detector()
        RUC._run_diffBragg_current(i_spot)
        RUC._set_background_plane(i_spot)
        RUC._extract_pixel_data()
        RUC._evaluate_averageI()

        img_minus = RUC.model_Lambda.copy()
        img_second_fdiff = (img_plus - 2*img0 + img_minus) / delta_h/delta_h
        all_delta_h2.append(delta_h**2)
        error2 = np.abs(d2 - img_second_fdiff).mean()
        all_error2.append( error2)

        if RUC.poisson_only:
            RUC._evaluate_log_averageI()
        else:
            RUC._evaluate_log_averageI_plus_sigma_readout()

        RUC._derivative_convenience_factors()
        Lminus = RUC._target_accumulate() #RUC._gaussian_target()

        finite_second_difference = (Lplus - 2*L + Lminus) / delta_h / delta_h
        Lerror2 = abs(finite_second_difference - dL_dtheta)
        all_Lerror2.append(Lerror2)

    if args.curvatures:
        print("Shift=%2.7g, image error=%2.7g,  L-error=%2.7g, image2=%2.7g, L2=%2.7g" %
              (delta_s, error, Lerror, error2, Lerror2))
    else:
        print("Shift=%2.7g, image error=%2.7g,  L-error=%2.7g" % (delta_s, error, Lerror))
    all_error.append(error)  # error in the image derivatives
    all_Lerror.append(Lerror)

from scipy.stats import linregress
l = linregress(all_delta_h, all_error)
assert l.slope > 0
assert l.rvalue > 0.999
assert l.pvalue < 1e-6

l = linregress(all_delta_h, all_Lerror)
assert l.slope > 0
assert l.rvalue > 0.999
assert l.pvalue < 1e-6

if args.curvatures:
    l = linregress(all_delta_h2, all_error2)
    assert l.slope > 0
    assert l.rvalue > 0.999
    assert l.pvalue < 1e-6

    l = linregress(all_delta_h2, all_Lerror2)
    assert l.slope > 0
    assert l.rvalue > 0.999
    assert l.pvalue < 1e-6

print("OK!")
#
#RUC.target_functional += RUC._target_accumulate()
## accumulate the per pixel derivatives
#RUC._background_derivatives(i_spot)
#RUC._Umatrix_derivatives()
#RUC._Bmatrix_derivatives()
#RUC._mosaic_parameter_m_derivatives()
#RUC._originZ_derivatives()
#RUC._spot_scale_derivatives()
#RUC._gain_factor_derivatives()
#RUC._Fcell_derivatives(i_spot)
## Done with derivative accumulation
#embed()
#
##for RUC._i_shot in RUC.shot_ids:
##    if RUC._i_shot in RUC.bad_shot_list:
##        continue
##    RUC.scale_fac = RUC._get_spot_scale(RUC._i_shot)
##
##    # TODO: Omatrix update? All crystal models here should have the same to_primitive operation, ideally
##    RUC._update_beams()
##    RUC._update_umatrix()
##    RUC._update_ucell()
##    RUC._update_ncells()
##    RUC._update_rotXYZ()
##    n_spots = len(RUC.NANOBRAGG_ROIS[RUC._i_shot])
##    for i_spot in range(n_spots):
##
##        RUC._panel_id = RUC.PANEL_IDS[RUC._i_shot][i_spot]
##
##        if RUC.verbose and i_spot % RUC.spot_print_stride == 0:
##            print("diffBragg: img %d/%d; spot %d/%d; panel %d" \
##                  % (RUC._i_shot + 1, RUC.n_shots, i_spot + 1, n_spots, RUC._panel_id), flush=True)
##
##        RUC.Imeas = RUC.ROI_IMGS[RUC._i_shot][i_spot]
##        RUC._update_dxtbx_detector()
##        RUC._run_diffBragg_current(i_spot)
##        RUC._set_background_plane(i_spot)
##        RUC._extract_pixel_data()
##        RUC._evaluate_averageI()
##
##        if RUC.poisson_only:
##            RUC._evaluate_log_averageI()
##        else:
##            RUC._evaluate_log_averageI_plus_sigma_readout()
##
##        RUC._derivative_convenience_factors()
##        RUC.target_functional += RUC._target_accumulate()
##
##        # accumulate the per pixel derivatives
##        RUC._background_derivatives(i_spot)
##        RUC._Umatrix_derivatives()
##        RUC._Bmatrix_derivatives()
##        RUC._mosaic_parameter_m_derivatives()
##        RUC._originZ_derivatives()
##        RUC._spot_scale_derivatives()
##        RUC._gain_factor_derivatives()
##        RUC._Fcell_derivatives(i_spot)
##        # Done with derivative accumulation
##
#RUC._mpi_aggregation()
#
#RUC._f = RUC.target_functional
#RUC._g = RUC.grad
#RUC.g = RUC.grad  # TODO why all these repeated definitions ?, RUC.g is needed by _verify_diag
#RUC._curvature_analysis()
#
## reset ROI pixels TODO: is this necessary
#RUC.D.raw_pixels *= 0
#
############################################
############################################
