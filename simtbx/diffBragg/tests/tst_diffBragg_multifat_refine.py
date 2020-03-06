from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
parser.add_argument("--detdist", action='store_true', help='perturb then refine the detdist')
parser.add_argument("--ncells", action='store_true', help='perturb then refine the ncells')
parser.add_argument("--bg", action='store_true', help='refine bg planes... ')
parser.add_argument("--fixscale", action='store_true', help='fix the scale')
parser.add_argument("--bmatrix", action='store_true')
parser.add_argument("--umatrix", action='store_true')
parser.add_argument("--nshots", default=1, type=int)
parser.add_argument("--curvatures", action='store_true')
parser.add_argument("--psf", action='store_true')
parser.add_argument("--gain", action='store_true')
args = parser.parse_args()

#if args.detdist:
#    raise NotImplementedError("Not yet refining detdist for single shots")

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
from simtbx.diffBragg.refiners.global_refiner import FatRefiner
from IPython import embed
from simtbx.diffBragg.refiners.crystal_systems import MonoclinicManager

# containers for FatRefine
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
ucell = (55, 65, 75, 90, 95, 90)
ucell2 = (55, 65, 75, 90, 95, 90)
if args.bmatrix:
    ucell2 = (55.1, 65.2, 74.9, 90, 94.9, 90)
symbol = "P121"

from simtbx.diffBragg.utils import  fcalc_from_pdb
miller_array = fcalc_from_pdb(resolution=2, wavelength=1, algorithm='fft', ucell=ucell, symbol=symbol)
Ncells_gt = 12, 12, 12

N_SHOTS = args.nshots

np.random.seed(3142019)
detdists_gt = np.random.normal(150,0.1, N_SHOTS)
offsets = np.random.uniform(1, 3, N_SHOTS) * np.random.choice([1,-1], N_SHOTS)
originZ_gt = {}
all_dets = []
for i_shot in range(N_SHOTS):

    # FIRST WE GENERATE SOME RANDOM IMAGES

    # generate a random raotation
    rotation = Rotation.random(num=1, random_state=100 + i_shot)[0]
    Q = rec(rotation.as_quat(), n=(4, 1))
    rot_ang, rot_axis = Q.unit_quaternion_as_axis_and_angle()

    # generate a small perturbation rotation
    perturb_rot_axis = np.random.random(3)
    perturb_rot_axis /= np.linalg.norm(perturb_rot_axis)
    perturb_rot_ang = 0
    if args.umatrix:
        perturb_rot_ang = np.random.choice([0.02, 0.03, 0.04, .05])  # degrees

    # make the ground truth crystal:
    a_real, b_real, c_real = sqr(uctbx.unit_cell(ucell).orthogonalization_matrix()).transpose().as_list_of_lists()
    C = Crystal(a_real, b_real, c_real, symbol)
    C.rotate_around_origin(rot_axis, rot_ang)

    # make the perturbed crystal model
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
    nbcryst.Ncells_abc = Ncells_gt  # ground truth Ncells

    nbcryst.miller_array = miller_array
    print("Ground truth ncells = %f" % (nbcryst.Ncells_abc[0]))

    # ground truth detector
    DET_gt = SimData.simple_detector(detdists_gt[i_shot], 0.177, (600, 600))
    originZ_gt[i_shot] = DET_gt[0].get_origin()[2]

    # initialize the simulator
    SIM = SimData()
    SIM.detector = DET_gt
    all_dets.append(DET_gt)

    # TODO get the detector model
    node = SIM.detector[0]
    node_d = node.to_dict()
    Origin = node_d["origin"][0], node_d["origin"][1], node_d["origin"][2]
    distance = Origin[2]
    print "Ground truth originZ=%f" % (SIM.detector[0].get_origin()[2])

    # TODO perturb the detector model
    # copy the detector and update the origin
    det2 = deepcopy(SIM.detector)
    # alter the detector distance by 2 mm
    #detz_offset = np.random.normal(0, 0.1)  # TODO make me random per shot (like GDVN would have)
    detz_offset = offsets[i_shot]
    node_d["origin"] = Origin[0], Origin[1], Origin[2] + detz_offset
    det2[0] = Panel.from_dict(node_d)

    shot_originZ_init[i_shot] = Origin[2]
    if args.detdist:
        shot_originZ_init[i_shot] = Origin[2]+detz_offset

    SIM.crystal = nbcryst
    SIM.instantiate_diffBragg(oversample=0)
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
    SIM.D.readout_noise_adu = 0
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

    print "Using oversample %d" % SIM.D.oversample

    # This is the ground truth image:
    img = SIM.D.raw_pixels.as_numpy_array()
    SIM.D.raw_pixels *= 0

    if args.psf:
        y = slice(450, 480,1)
        x = slice(650, 670, 1)
        print "PSF max discrepancy: %f" % abs(img_pre_psf[y,x]- img[y,x]).max()

    # Simulate the perturbed image for comparison
    # perturbed detector:
    if args.detdist:
        SIM.detector = det2
        SIM.D.update_dxtbx_geoms(det2, SIM.beam.nanoBragg_constructor_beam, 0)
        print ("Modified originZ=%f" % (det2[0].get_origin()[2]))
    # perturbed crystal:
    SIM.D.Bmatrix = C2.get_B()
    SIM.D.Umatrix = C2.get_U()
    nbcryst.dxtbx_crystal = C2
    if args.ncells:
        Ncells_abc2 = 14, 14, 14
        nbcryst.Ncells_abc = Ncells_abc2
        SIM.D.set_value(9, Ncells_abc2[0])
        print ("Modified Ncells=%f" % Ncells_abc2[0])

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
    spot_roi, tilt_abc = utils.process_simdata(spots, img, thresh=20, plot=args.plot, shoebox_sz=30)

    UcellMan = MonoclinicManager(
        a=ucell2[0],
        b=ucell2[1],
        c=ucell2[2],
        beta=ucell2[4]*np.pi/180.)

    if args.gain:
        img = img*1.1

    # TODO: the following need to be embedded in the refiner init function..
    nspot = len(spot_roi)
    nspot_per_shot[i_shot] = nspot

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
        SIM.D.free_all()  # CLEANGIUAGE

if args.detdist:
    SIM.D.oversample_omega = False  # necessary to refine detector distance

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

nrotation_param = 3*N_SHOTS
nscale_param = 1*N_SHOTS
ntilt_param = 0
for i_shot in range(N_SHOTS):
    ntilt_param += 3 * nspot_per_shot[i_shot]

ndetz_param = len(detdists_gt)
n_local_unknowns = nrotation_param + nscale_param + ntilt_param + ndetz_param

nucell_param = len(UcellMan.variables)
n_ncell_param = 1
nfcell_param = len(idx_from_asu)
ngain_param = 1

n_global_unknowns = nucell_param + nfcell_param + ngain_param + n_ncell_param
n_total_unknowns = n_local_unknowns + n_global_unknowns

RUC = FatRefiner(
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
    global_ncells=True,
    global_ucell=True,
    global_originZ=False,
    shot_originZ_init=shot_originZ_init,
    sgsymbol=symbol)

RUC.idx_from_asu = idx_from_asu
RUC.asu_from_idx = asu_from_idx
RUC.refine_background_planes = args.bg
RUC.refine_Umatrix = args.umatrix
RUC.refine_Bmatrix = args.bmatrix
RUC.refine_ncells = args.ncells
RUC.refine_crystal_scale = not args.fixscale
RUC.refine_Fcell = False
RUC.refine_detdist = args.detdist
RUC.refine_gain_fac = args.gain

RUC.max_calls = 3000
RUC.trad_conv_eps = 1e-2
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
RUC.has_pre_cached_roi_data = True
RUC.S.D.update_oversample_during_refinement = False
RUC.use_curvatures = False
RUC.use_curvatures_threshold = 4
RUC.bg_offset_positive = args.bg
RUC.bg_offset_only = args.bg
RUC.calc_curvatures = args.curvatures
RUC.poisson_only = True
RUC.verbose = True
RUC.big_dump = True
RUC.gt_ncells = Ncells_gt[0]
RUC.originZ_gt = originZ_gt
RUC.gt_ucell = ucell[0], ucell[1], ucell[2], ucell[4]
RUC.testing_mode = True
RUC.run(setup_only=False)
RUC.run(setup_only=True)
if RUC.hit_break_to_use_curvatures:
    RUC.num_positive_curvatures = 0
    RUC.use_curvatures = True
    RUC.run(setup=False)



#RUC.calc_func = True
#RUC.compute_functional_and_gradients()
#
#def func(x, RUC):
#    print("F: det dist %f" % RUC.x[-3])
#    RUC.calc_func = True
#    RUC.x = flex.double(x)
#    f, g = RUC.compute_functional_and_gradients()
#    return f
#
#
#def fprime(x, RUC):
#    print("G: det dist %f" % RUC.x[-3])
#    RUC.calc_func = False
#    RUC.x = flex.double(x)
#    RUC.x = flex.double(x)
#    f, g = RUC.compute_functional_and_gradients()
#    return 1*g.as_numpy_array()
#
#
#from scipy.optimize import fmin_l_bfgs_b
#
#bounds = [(-np.inf, np.inf)]*RUC.n
#bounds[-11] = -.1*np.pi/180, .1*np.pi/180  # rotX
#bounds[-10] = -.1*np.pi/180, .1*np.pi/180  # roty
#bounds[-9] = -.1*np.pi/180, .1*np.pi/180  # rotZ
#bounds[-8] = 50, 60  # a
#bounds[-7] = 60, 70  #  b
#bounds[-6] = 70, 80  # c
#bounds[-5] = 93*np.pi/180, 97*np.pi/180.  # beta
#bounds[-4] = 7, 30  # ncells
#bounds[-3] = -170, -150  # detdist
#bounds[-2] = 1, 1  # gain
#bounds[-1] = 1, 1  # scale
#
#print("GO!")
##out = fmin_l_bfgs_b(func=func, x0=np.array(RUC.x),
##                    fprime=fprime,args=[RUC]) #, bounds=bounds)
#out = fmin_l_bfgs_b(func=func, factr=1000, x0=np.array(RUC.x),
#                    fprime=fprime, maxls=100,
#                    pgtol=1e-7,
#                    args=(RUC,),
#                    bounds=bounds)

