
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
parser.add_argument("--refine", action='store_true')
parser.add_argument("--crystalsystem", default='tetragonal',
                    choices=["monoclinic", "tetragonal"])
args = parser.parse_args()
is_tet = args.crystalsystem == "tetragonal"
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.stats import pearsonr
import pylab as plt

from cctbx import sgtbx
from rstbx.symmetry.constraints import parameter_reduction
from scitbx.matrix import sqr
from scitbx.matrix import rec
from simtbx.diffBragg.nanoBragg_crystal import nanoBragg_crystal
from cctbx import uctbx
from simtbx.diffBragg.sim_data2 import SimData
from simtbx.diffBragg import utils
from dxtbx.model.crystal import Crystal
from IPython import embed

#TODO parameters shifts should be percentages

# STEP 1:
# make a crystal and orient it randomly

if is_tet:
    ucell = (55, 55, 77, 90, 90, 90)
    symbol = "P43212"
else:  # args.crystalsystem == "monoclinic"
    ucell = (55, 65, 75, 90, 95, 90)
    symbol = "P121"

a_real, b_real, c_real = sqr(uctbx.unit_cell(ucell).orthogonalization_matrix()).transpose().as_list_of_lists()
#a_real, b_real, c_real = np.reshape(uctbx.unit_cell(ucell).orthogonalization_matrix(), (3,3)).T   #transpose().as_list_of_lists()
C = Crystal(a_real, b_real, c_real, symbol)

# random raotation
rotation = Rotation.random(num=1, random_state=101)[0]
Q = rec(rotation.as_quat(), n=(4, 1))
rot_ang, rot_axis = Q.unit_quaternion_as_axis_and_angle()
C.rotate_around_origin(rot_axis, rot_ang)
a_rot, b_rot, c_rot = C.get_real_space_vectors()

# make the parameter reduction manager
point_group = sgtbx.space_group_info(symbol=symbol).group().build_derived_point_group()
S = parameter_reduction.symmetrize_reduce_enlarge(point_group)
S.set_orientation(np.hstack((a_real, b_real, c_real)), length_unit=1e-10)
X = S.forward_independent_parameters()
dX = S.forward_gradients()  # gradients in meters
n_ucell_params = len(X)


# STEP3: sanity check
B = S.backward_orientation(independent=X)
print("Number of parameters = %d" % n_ucell_params)
print("Test congruency between parameter reduction module and dxtbx...")
assert np.allclose(C.get_B(), B.reciprocal_matrix())
#Bstar = B.reciprocal_matrix()
print("OK!")

if is_tet:
    UcellMan = utils.TetragonalManager(a=ucell[0], c=ucell[2])
else:
    UcellMan = utils.MonoclinicManager(a=ucell[0], b=ucell[1],
                                       c=ucell[2], beta=ucell[4]*np.pi/180)

assert np.allclose(UcellMan.B_recipspace, C.get_B())

# STEP4:
# make a nanoBragg crystal to pass to diffBragg
nbcryst = nanoBragg_crystal()
nbcryst.dxtbx_crystal = C
nbcryst.n_mos_domains = 1
nbcryst.thick_mm = 0.01
nbcryst.Ncells_abc = (7, 7, 7)

# STEP5: make an instance of diffBRagg, use the simData wrapper
SIM = SimData()
# overwrite the default detector with a smaller pixels one
SIM.detector = SimData.simple_detector(298, 0.1, (700, 700))
# FIXME: determine why setting detdist > 298 causes oversample to diverge...
SIM.crystal = nbcryst
SIM.instantiate_diffBragg(oversample=0)
# D is an instance of diffBragg with reasonable parameters
# and our dxtbx crystal created above
D = SIM.D

# STEP6:
# initialize the derivative managers for the unit cell parameters
for i_param in range(n_ucell_params):
    D.refine(3+i_param)
D.initialize_managers()
for i in range(n_ucell_params):
    D.set_ucell_derivative_matrix(3+i, UcellMan.derivative_matrices[i])
D.initialize_managers()

roi = ((0, 699), (0, 699))
rX = slice(roi[0][0], roi[0][1], 1)
rY = slice(roi[1][0], roi[1][1], 1)
D.region_of_interest = roi

# STEP7:
# compute the scattering and its derivative
D.add_diffBragg_spots()
img = D.raw_pixels_roi.as_numpy_array()
# reset all pixel values
D.raw_pixels *= 0
D.raw_pixels_roi *= 0

derivs = []
for i_param in range(n_ucell_params):
    analy_deriv = SIM.D.get_derivative_pixels(3+i_param).as_numpy_array()
    derivs.append(analy_deriv)

# STEP8
# iterate over the parameters and do a finite difference test for each one
# parameter shifts:
shifts = [1e-3 * (2**i) for i in range(1, 12, 2)]

if not args.refine:
    for i_param in range(n_ucell_params):
        #if i_param < 3:
        #    continue
        analy_deriv = derivs[i_param]
        diffs = []
        for i_shift, param_shift in enumerate(shifts):

            #uc2 = C2.get_unit_cell().parameters()
            #Tet.a = a2_real[0]
            #Tet.c = c2_real[2]
            if is_tet:
                var = [ucell[0], ucell[2]]
            else:
                var = [ucell[0], ucell[1], ucell[2], ucell[4]*np.pi/180]

            name = UcellMan.variable_names[i_param]
            if "alpha" in name or "beta" in name or "gamma" in name:
                param_shift = param_shift * 1e-4  # rescale the radian parameters

            var[i_param] += param_shift
            UcellMan.variables = var

            D.Bmatrix = UcellMan.B_recipspace
            D.add_diffBragg_spots()

            img2 = D.raw_pixels_roi.as_numpy_array()

            # reset for next computation
            D.raw_pixels_roi *= 0
            D.raw_pixels *= 0

            finite_deriv = (img2-img) / param_shift

            y = slice(340, 360)
            x = slice(240, 260)
            bragg = img > 0.5

            ave_error = np.abs(finite_deriv[bragg] - analy_deriv[bragg]).mean()
            print "\tAverage error=%f; parameter shift h=%f" % (ave_error, abs(param_shift))

            r = pearsonr(analy_deriv[bragg].ravel(), finite_deriv[bragg].ravel())[0]
            diffs.append(r)
            if args.plot:
                plt.subplot(121)
                plt.imshow(finite_deriv)
                plt.title("finite diff")
                plt.subplot(122)
                plt.imshow(analy_deriv)
                plt.title("analytical")
                plt.draw()
                plt.suptitle("Shift %d / %d"
                             % (i_shift+1, len(shifts)))
                plt.pause(0.8)

        if args.plot:
            plt.close()
            plt.plot(shifts, diffs, 'o')
            title = "Unit cell parameter %d / %d" % (i_param+1, n_ucell_params)
            plt.title(title + "\nPearson corr between finite diff and analytical")
            plt.xlabel("unit cell shifts")
            plt.ylabel("Pearson corr")
            plt.show()

        # verify a high correlation for the smallest parameter shift
        print("Check high pearson R between analytical and finite diff")
        print("Pearson correlection at smallest parameter shift=%f" % diffs[0])
        assert(diffs[0] > .98), "%f" % diffs[0]
        # check monotonic decrease
        print("Fit polynomial and check monotonic decrease")
        trend = np.polyval(np.polyfit(shifts, diffs, 2), shifts)
        assert np.all(np.diff(zip(trend[:-1], trend[1:]), axis=1) <= 0)
    print("OK!")


if is_tet:
    ucell2 = (55.5, 55.5, 76.1, 90, 90, 90)
else:  # args.crystalsystem == "monoclinic"
    ucell2 = (55.2, 66, 74, 90, 94.3, 90)

a2_real, b2_real, c2_real = sqr(uctbx.unit_cell(ucell2).orthogonalization_matrix()).transpose().as_list_of_lists()
#a_real, b_real, c_real = np.reshape(uctbx.unit_cell(ucell).orthogonalization_matrix(), (3,3)).T   #transpose().as_list_of_lists()
C2 = Crystal(a2_real, b2_real, c2_real, symbol)

print ("Ensure the U matrices are the same for the ground truth and perturbation")
np.random.seed(1)
axis = np.random.random(3)
axis /= np.linalg.norm(axis)
C2.rotate_around_origin(rot_axis, rot_ang)
assert np.allclose(C.get_U(), C2.get_U())


# Setup the simulation and create a realistic image
# with background and noise
# <><><><><><><><><><><><><><><><><><><><><><><><><>
nbcryst.dxtbx_crystal = C   # simulate ground truth
nbcryst.thick_mm = 0.5
nbcryst.Ncells_abc = 12, 12, 12

SIM = SimData()
SIM.detector = SimData.simple_detector(150, 0.1, (512, 512))
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
# this is the ground truth image:
img = SIM.D.raw_pixels.as_numpy_array()
SIM.D.raw_pixels *= 0

# Simulate the perturbed image for comparison
SIM.D.Bmatrix = C2.get_B()
SIM.D.add_diffBragg_spots()
SIM._add_background()
SIM._add_noise()
img_pet = SIM.D.raw_pixels.as_numpy_array()
SIM.D.raw_pixels *= 0

# NOTE NEED TO DO SPOT FINDING AND WHAT NOT
# spot_rois, abc_init, img, SimData_instance,
# locate the strong spots and fit background planes
# <><><><><><><><><><><><><><><><><><><><><><><><><>
spot_data = utils.get_spot_data(spots, thresh=19)
ss_spot, fs_spot = map(np.array, zip(*spot_data["maxIpos"]))  # slow/fast  scan coords of strong spots
num_spots = len(ss_spot)
if args.plot:
    plt.imshow(img, vmax=200)
    plt.plot(fs_spot, ss_spot, 'o', mfc='none', mec='r')
    plt.title("Simulated image with strong spots marked")
    plt.show()

is_bg_pixel = np.ones(img.shape, bool)
for bb_ss, bb_fs in spot_data["bboxes"]:
    is_bg_pixel[bb_ss, bb_fs] = False

# now fit tilting planes
shoebox_sz = 20
tilt_abc = np.zeros((num_spots, 3))
spot_roi = np.zeros((num_spots, 4), int)
if args.plot:
    patches = []
img_shape = SIM.detector[0].get_image_size()
for i_spot, (x_com, y_com) in enumerate(zip(fs_spot, ss_spot)):
    i1 = int(max(x_com - shoebox_sz / 2., 0))
    i2 = int(min(x_com + shoebox_sz / 2., img_shape[0]-1))
    j1 = int(max(y_com - shoebox_sz / 2., 0))
    j2 = int(min(y_com + shoebox_sz / 2., img_shape[1]-1))

    shoebox_img = img[j1:j2, i1:i2]
    shoebox_mask = is_bg_pixel[j1:j2, i1:i2]

    tilt, bgmask, coeff = utils.tilting_plane(
        shoebox_img,
        mask=shoebox_mask,  # mask specifies which spots are bg pixels...
        zscore=2)

    tilt_abc[i_spot] = coeff[1], coeff[2], coeff[0]  # store as fast-scan coeff, slow-scan coeff, offset coeff

    spot_roi[i_spot] = i1, i2, j1, j2
    if args.plot:
        R = plt.Rectangle(xy=(x_com-shoebox_sz/2, y_com-shoebox_sz/2.),
                      width=shoebox_sz,
                      height=shoebox_sz,
                      fc='none', ec='r')
        patches.append(R)

if args.plot:
    patch_coll = plt.mpl.collections.PatchCollection(patches, match_original=True)
    plt.imshow(img, vmin=0, vmax=200)
    plt.gca().add_collection(patch_coll)
    plt.show()

from simtbx.diffBragg.refiners import RefineUcell

if is_tet:
    UcellMan.variables = [ucell2[0], ucell2[2]]
else:
    UcellMan.variables = [ucell2[0], ucell2[1], ucell2[2], ucell2[4]*np.pi/180]

nbcryst.dxtbx_crystal = C2
SIM.crystal = nbcryst

RUC = RefineUcell(
    spot_rois=spot_roi,
    abc_init=tilt_abc,
    img=img,
    SimData_instance=SIM,
    plot_images=args.plot,
    ucell_manager=UcellMan)

RUC.trad_conv = True
RUC.trad_conv_eps = 1e-5
RUC.run()

