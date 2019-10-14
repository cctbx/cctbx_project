from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
parser.add_argument("--detdist", action='store_true', help='perturb then refine the detdist')
parser.add_argument("--ncells", action='store_true', help='perturb then refine the ncells')
parser.add_argument("--bmatrix", action='store_true')
parser.add_argument("--umatrix", action='store_true')
parser.add_argument("--curvatures", action='store_true')
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
    perturb_rot_ang = .15  # degrees

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
Ncells_gt = 12,12,12
nbcryst.Ncells_abc = Ncells_gt  # ground truth Ncells
print("Ground truth ncells = %f" % (nbcryst.Ncells_abc[0]))

SIM = SimData()
SIM.detector = SimData.simple_detector(150, 0.1, (512, 512))

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
detz_offset = 0.3
node_d["origin"] = Origin[0], Origin[1], Origin[2] + detz_offset
det2[0] = Panel.from_dict(node_d)

SIM.crystal = nbcryst
SIM.instantiate_diffBragg(oversample=0)
SIM.D.nopolar = False
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
print SIM.D.oversample

# This is the ground truth image:
img = SIM.D.raw_pixels.as_numpy_array()
SIM.D.raw_pixels *= 0

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

RUC = RefineAll(
    spot_rois=spot_roi,
    abc_init=tilt_abc,
    img=1.1*img,
    SimData_instance=SIM,
    plot_images=args.plot,
    plot_residuals=True,
    ucell_manager=UcellMan)

RUC.trad_conv = True
RUC.refine_detdist = args.detdist
RUC.refine_background_planes = False
RUC.refine_Umatrix = args.umatrix
RUC.refine_Bmatrix = args.bmatrix
RUC.refine_ncells = args.ncells
RUC.use_curvatures = args.curvatures
RUC.refine_crystal_scale = True
RUC.refine_gain_fac = True
RUC.plot_stride = 10
RUC.plot_residuals = args.plot
RUC.trad_conv_eps = 1e-5
RUC.max_calls = 300
#RUC._setup()
#RUC._cache_roi_arrays()
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
    if ang > 0 and perturb_rot_ang > 0:
        assert np.linalg.norm(np.round(ax, 2) + np.round(perturb_rot_axis, 2)) < 0.075
        assert abs(ang - perturb_rot_ang) < 0.01
    assert final_Umat_norm < 1e-1*init_Umat_norm

if args.ncells:
    assert round(RUC.x[-4]) == Ncells_gt[0]

print("OK")

