from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
parser.add_argument("--detdist", action='store_true', help='perturb then refine the detdist')
parser.add_argument("--ncells", action='store_true', help='perturb then refine the ncells')
parser.add_argument("--bmatrix", action='store_true')
parser.add_argument("--umatrix", action='store_true')
parser.add_argument("--fixscale", action='store_true', help='fix the scale')
parser.add_argument("--curvatures", action='store_true')
parser.add_argument("--psf", action='store_true')
parser.add_argument("--gain", action='store_true')
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
from scitbx.matrix import sqr

from simtbx.nanoBragg import shapetype
from simtbx.diffBragg.utils import fcalc_from_pdb
from simtbx.diffBragg.nanoBragg_crystal import nanoBragg_crystal
from simtbx.diffBragg.sim_data import SimData
from simtbx.diffBragg import utils
from simtbx.diffBragg.refiners import RefineAll
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
nbcryst = nanoBragg_crystal()
nbcryst.dxtbx_crystal = C   # simulate ground truth
nbcryst.thick_mm = 0.1
Ncells_gt = 12, 12, 12
nbcryst.Ncells_abc = Ncells_gt  # ground truth Ncells
nbcryst.miller_array = fcalc_from_pdb(resolution=2, wavelength=1,algorithm='fft', ucell=ucell, symbol=symbol)
print("Ground truth ncells = %f" % (nbcryst.Ncells_abc[0]))

SIM = SimData()
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

print ("Using oversample %d" % SIM.D.oversample)

# This is the ground truth image:
img = SIM.D.raw_pixels.as_numpy_array()
SIM.D.raw_pixels *= 0

if args.psf:
    y = slice(450, 480,1)
    x = slice(650, 670, 1)
    print ("PSF max discrepancy: %f" % abs(img_pre_psf[y,x]- img[y,x]).max())

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

#print("I got %s spots to process!" % spot_roi.shape[0])
#n_kept = 30
#np.random.seed(2)
#idx = np.random.permutation(spot_roi.shape[0])[:n_kept]
#spot_roi = spot_roi[idx]
#tilt_abc = tilt_abc[idx]
#print ("I kept %d spots!" % tilt_abc.shape[0])


UcellMan = MonoclinicManager(
    a=ucell2[0],
    b=ucell2[1],
    c=ucell2[2],
    beta=ucell2[4]*np.pi/180.)

init_Umat_norm = np.abs(np.array(C2.get_U()) - np.array(C.get_U())).sum()
init_Bmat_norm = np.abs(np.array(C2.get_B()) - np.array(C.get_B())).sum()

if args.gain:
    img = img*1.1

RUC = RefineAll(
    spot_rois=spot_roi,
    abc_init=tilt_abc,
    img=img,
    SimData_instance=SIM,
    plot_images=args.plot,
    plot_residuals=True,
    ucell_manager=UcellMan)


#from simtbx.diffBragg.refiners.global_refiner import FatRefiner
#
#
## TODO: the following need to be embedded in the refiner init function..
#nspot = len(spot_roi)
#
#
#
#
#nanoBragg_rois = []  # special nanoBragg format
#xrel, yrel, roi_imgs = [], [], []
#xcom, ycom = [],[]
#for i_roi, (x1, x2, y1, y2) in enumerate(spot_roi):
#    nanoBragg_rois.append(((x1, x2), (y1, y2)))
#    yr, xr = np.indices((y2 - y1 + 1, x2 - x1 + 1))
#    xrel.append(xr)
#    yrel.append(yr)
#    roi_imgs.append(img[y1:y2 + 1, x1:x2 + 1])
#
#    xcom.append(.5*(x1 + x2))
#    ycom.append(.5*(x1 + x2))
#
#q_spot = utils.x_y_to_q(xcom, ycom, SIM.detector, SIM.beam.nanoBragg_constructor_beam)
#Ai = sqr(SIM.crystal.dxtbx_crystal.get_A()).inverse()
#Ai = Ai.as_numpy_array()
#HKL = np.dot(Ai, q_spot.T)
#HKLi = [np.ceil(h - 0.5).astype(int) for h in HKL]
#HKLi = [tuple(x) for x in np.vstack(HKLi).T]
#Hi_asu = utils.map_hkl_list(HKLi, anomalous_flag=True, symbol=symbol)
#
#nrotation_param = 3
#nscale_param = 1
#ntilt_param = 3*nspot
#n_local_unknowns = nrotation_param + nscale_param + ntilt_param
#
#nucell_param = len(UcellMan.variables)
#n_ncell_param = 1
#nfcell_param = len(Hi_asu)
#ngain_param = 1
#ndetz_param = 1
#
#n_global_unknowns = nucell_param + nfcell_param + ngain_param + ndetz_param + n_ncell_param
#n_total_unknowns = n_local_unknowns + n_global_unknowns
#
#
#RUC = FatRefiner(
#    n_total_params=n_total_unknowns,
#    n_local_params=n_local_unknowns,
#    n_global_params=n_global_unknowns,
#    local_idx_start=0,
#    shot_ucell_managers={0: UcellMan},
#    shot_rois={0: spot_roi},
#    shot_nanoBragg_rois={0: nanoBragg_rois},
#    shot_roi_imgs={0: roi_imgs},
#    shot_spectra={0: SIM.beam.spectrum},
#    shot_crystal_GTs={0: C},
#    shot_crystal_models={0: SIM.crystal.dxtbx_crystal},
#    shot_xrel={0: xrel},
#    shot_yrel={0: yrel},
#    shot_abc_inits={0: tilt_abc},
#    shot_asu={0: Hi_asu},
#    global_param_idx_start=n_local_unknowns,
#    shot_panel_ids={0: [0]*nspot},
#    log_of_init_crystal_scales=None,
#    all_crystal_scales=None,
#    perturb_fcell=False,
#    global_ncells=True,
#    global_ucell=True,
#    sgsymbol=symbol)
#
##TODO make this part of class init:
#idx_from_asu = {h: i for i, h in enumerate(set(Hi_asu))}
#asu_from_idx = {i: h for i, h in enumerate(set(Hi_asu))}
#
#RUC.idx_from_asu = idx_from_asu
#RUC.asu_from_idx = asu_from_idx
#
#RUC.refine_background_planes = False
#RUC.refine_Umatrix = args.umatrix
#RUC.refine_Bmatrix = args.bmatrix
#RUC.refine_ncells = args.ncells
#RUC.refine_crystal_scale = False
#RUC.refine_crystal_scale = True
#RUC.refine_Fcell = False
#RUC.refine_detdist = args.detdist
#RUC.refine_gain_fac = args.gain
#
#RUC.max_calls = 3000
#RUC.trad_conv_eps = 1e-2
#RUC.trad_conv = True
#RUC.trial_id = 0
#
#
#RUC.plot_stride = 10
#RUC.plot_residuals = False
#RUC.plot_images = args.plot
#RUC.setup_plots()
#
#RUC.refine_rotZ = True
#RUC.request_diag_once = False
#RUC.S = SIM
#RUC.has_pre_cached_roi_data = True
#RUC.S.D.update_oversample_during_refinement = False
#RUC.use_curvatures = False  # args.curvatures
#RUC.use_curvatures_threshold = 1
#RUC.calc_curvatures = args.curvatures
#RUC.poisson_only = True
#RUC.verbose = True
#RUC.big_dump = True

#============




RUC.trad_conv = True
RUC.refine_detdist = args.detdist
RUC.refine_background_planes = False
RUC.refine_Umatrix = args.umatrix
RUC.refine_Bmatrix = args.bmatrix
RUC.refine_ncells = args.ncells
RUC.use_curvatures = args.curvatures
RUC.refine_crystal_scale =  not args.fixscale
RUC.refine_gain_fac = False
RUC.plot_stride = 10
RUC.plot_residuals = args.plot
RUC.trad_conv_eps = 1e-5
RUC.max_calls = 3000
#RUC._setup()
#RUC._cache_roi_arrays()
RUC.testing_mode = True
RUC.gt_ucell = ucell[0], ucell[1], ucell[2], ucell[4]
RUC.CRYSTAL_MODELS = {0: RUC.S.crystal.dxtbx_crystal}
RUC.CRYSTAL_GT = {0: C}
RUC.symbol = symbol
RUC.gt_ncells = Ncells_gt[0]
RUC.run()



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

ang, ax = RUC.get_correction_misset(as_axis_angle_deg=True)
if ang > 0:
    C2.rotate_around_origin(ax, ang)
C2.set_B(RUC.get_refined_Bmatrix())

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

if args.bmatrix:
    assert err_ref < 1e-1 * err_init

# NOTE, this test might change, e.g. angle could be negative and axis could be the same...
if args.umatrix:
    # the initial perturbation matrix:
    R1 = rec(perturb_rot_axis, (3, 1)).axis_and_angle_as_r3_rotation_matrix(perturb_rot_ang, deg=True)
    # restoring purturbation applied after refinement:
    R2 = ax.axis_and_angle_as_r3_rotation_matrix(ang, deg=True)

    # the hope is that the refined R2 cancels the effect of R1
    # hence, the product R1 and R2 should be ~ Identity
    I = np.reshape(R1 * R2, (3, 3))
    assert np.all(np.round(I - np.eye(3), 3) == np.zeros((3, 3)))
    assert final_Umat_norm < 1e-1*init_Umat_norm

if args.ncells:
    ncells_val = np.exp(RUC.x[RUC.ncells_xpos[0]]) + 3
    #ncells_val = RUC.x[RUC.ncells_xpos[0]]
    assert np.round(ncells_val) == Ncells_gt[0]

print("OK")

