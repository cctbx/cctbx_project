
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
parser.add_argument("--refine", action='store_true')
parser.add_argument("--crystalsystem", default='tetragonal',
                    choices=["monoclinic", "tetragonal"])
parser.add_argument("--curvatures", action='store_true')
args = parser.parse_args()
is_tet = args.crystalsystem == "tetragonal"
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.stats import pearsonr
import pylab as plt

from scitbx.matrix import sqr
from scitbx.matrix import rec
from simtbx.diffBragg.nanoBragg_crystal import nanoBragg_crystal
from cctbx import uctbx
from simtbx.diffBragg.sim_data import SimData
from simtbx.diffBragg import utils
from simtbx.diffBragg.refiners import crystal_systems
from dxtbx.model.crystal import Crystal
from simtbx.diffBragg.refiners import RefineUcell

#TODO parameters shifts should be percentages

# STEP 1:
# make a crystal and orient it randomly

if is_tet:
    ucell = (55, 55, 77, 90, 90, 90)
    symbol = "P43212"
else:  # args.crystalsystem == "monoclinic"
    ucell = (70, 60, 50, 90.0, 110, 90.0)
    symbol = "C121"

a_real, b_real, c_real = sqr(uctbx.unit_cell(ucell).orthogonalization_matrix()).transpose().as_list_of_lists()
C = Crystal(a_real, b_real, c_real, symbol)

# random raotation
rotation = Rotation.random(num=1, random_state=101)[0]
Q = rec(rotation.as_quat(), n=(4, 1))
rot_ang, rot_axis = Q.unit_quaternion_as_axis_and_angle()
C.rotate_around_origin(rot_axis, rot_ang)
a_rot, b_rot, c_rot = C.get_real_space_vectors()

# STEP3:
# create the unit cell parameter manager
if is_tet:
    UcellMan = crystal_systems.TetragonalManager(a=ucell[0], c=ucell[2])
else:
    UcellMan = crystal_systems.MonoclinicManager(a=ucell[0], b=ucell[1],
                                       c=ucell[2], beta=ucell[4]*np.pi/180)

n_ucell_params = len(UcellMan.variables)

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
SIM.detector = SimData.simple_detector(300, 0.1, (700, 700))
SIM.crystal = nbcryst
SIM.instantiate_diffBragg(oversample=0, verbose=0)
# D is an instance of diffBragg with reasonable parameters
# and our dxtbx crystal created above
D = SIM.D
D.progress_meter = True

# STEP6:
# initialize the derivative managers for the unit cell parameters
for i_param in range(n_ucell_params):
    D.refine(3+i_param)
D.initialize_managers()
for i in range(n_ucell_params):
    D.set_ucell_derivative_matrix(3+i, UcellMan.derivative_matrices[i])
    if args.curvatures:
        D.set_ucell_second_derivative_matrix(3 + i, UcellMan.second_derivative_matrices[i])
D.initialize_managers()

roi = ((0, 699), (0, 699))
rX = slice(roi[0][0], roi[0][1], 1)
rY = slice(roi[1][0], roi[1][1], 1)
D.region_of_interest = roi

# STEP7:
# compute the scattering and its derivative
print("Adding diffBragg spots")
D.add_diffBragg_spots()
print("Done!")
img = D.raw_pixels_roi.as_numpy_array()
# reset all pixel values
D.raw_pixels *= 0
D.raw_pixels_roi *= 0

derivs = []
second_derivs = []
for i_param in range(n_ucell_params):
    analy_deriv = D.get_derivative_pixels(3+i_param).as_numpy_array()
    derivs.append(analy_deriv)
    if args.curvatures:
        second_derivs.append(D.get_second_derivative_pixels(3+i_param).as_numpy_array())

# STEP8
# iterate over the parameters and do a finite difference test for each one
# parameter shifts:
shifts = [1e-3 * (2**i) for i in range(1, 12, 2)]

if not args.refine:
    for i_param in range(n_ucell_params):
        analy_deriv = derivs[i_param]
        diffs = []
        diffs2 = []
        for i_shift, param_shift in enumerate(shifts):

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

            img_forward = D.raw_pixels_roi.as_numpy_array()

            # reset for next computation
            D.raw_pixels_roi *= 0
            D.raw_pixels *= 0

            if args.curvatures:
                # estimate the second derivative
                var[i_param] = var[i_param] - 2*param_shift  # do the backwards finite deriv
                UcellMan.variables = var

                D.Bmatrix = UcellMan.B_recipspace
                D.add_diffBragg_spots()

                img_backward = D.raw_pixels_roi.as_numpy_array()

                # reset for next computation
                D.raw_pixels_roi *= 0
                D.raw_pixels *= 0

            finite_deriv = (img_forward-img) / param_shift

            if second_derivs:
                finite_second_deriv = (img_forward - 2*img + img_backward) / param_shift / param_shift

            bragg = img > 0.5

            ave_error = np.abs(finite_deriv[bragg] - analy_deriv[bragg]).mean()

            r = pearsonr(analy_deriv[bragg].ravel(), finite_deriv[bragg].ravel())[0]
            diffs.append(r)

            print ("\tAverage error=%f; parameter shift h=%f" % (ave_error, abs(param_shift)))
            if args.curvatures:
                ave_error2 = np.abs(finite_second_deriv[bragg] - second_derivs[i_param][bragg]).mean()
                print("\tsecond derivative Average error=%f; parameter shift squared h^2=%f"
                      % (ave_error2, abs(param_shift)**2))
                r2 = pearsonr(second_derivs[i_param][bragg].ravel(), finite_second_deriv[bragg].ravel())[0]
                diffs2.append(r2)
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
                if args.curvatures:
                    plt.subplot(121)
                    plt.imshow(finite_second_deriv)
                    plt.title("finite second diff")
                    plt.subplot(122)
                    plt.imshow(second_derivs[i_param])
                    plt.title("analytical")
                    plt.draw()
                    plt.suptitle("Shift %d / %d"
                                 % (i_shift + 1, len(shifts)))
                    plt.pause(0.8)

        if args.plot:
            plt.close()
            plt.plot(shifts, diffs, 'o')
            title = "Unit cell parameter %d / %d" % (i_param+1, n_ucell_params)
            plt.title(title + "\nPearson corr between finite deriv and analytical")
            plt.xlabel("unit cell shifts")
            plt.ylabel("Pearson corr")
            plt.show()
            if args.curvatures:
                plt.close()
                plt.plot(np.array(shifts)**2, diffs2, 'o')
                title = "Unit cell parameter %d / %d" % (i_param + 1, n_ucell_params)
                plt.title(title + "\nPearson corr between finite second deriv and analytical")
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
        if args.curvatures:
            assert (diffs2[0] > .99), "%f" % diffs2[0]
    print("OK!")

if is_tet:
    ucell2 = (55.5, 55.5, 76.1, 90, 90, 90)
else:  # args.crystalsystem == "monoclinic"
    ucell2 = (71, 59, 51, 90.0, 111, 90.0)

a2_real, b2_real, c2_real = sqr(uctbx.unit_cell(ucell2).orthogonalization_matrix()).transpose().as_list_of_lists()
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
SIM.detector = SimData.simple_detector(150, 0.1, (1024, 1024))
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
spot_roi, tilt_abc = utils.process_simdata(spots, img, thresh=20, plot=args.plot)

if is_tet:
    UcellMan.variables = [ucell2[0], ucell2[2]]
else:
    UcellMan.variables = [ucell2[0], ucell2[1], ucell2[2], ucell2[4]*np.pi/180]

nbcryst.dxtbx_crystal = C2
SIM.crystal = nbcryst
SIM.D.Omatrix = nbcryst.Omatrix
SIM.D.Bmatrix = C2.get_B()
SIM.D.Umatrix = C2.get_U()

RUC = RefineUcell(
    spot_rois=spot_roi,
    abc_init=tilt_abc,
    img=img,
    SimData_instance=SIM,
    plot_images=args.plot,
    ucell_manager=UcellMan)

RUC.trad_conv = True
RUC.trad_conv_eps = 1e-4
RUC.max_calls = 150
RUC.use_curvatures = args.curvatures
RUC.run()

if is_tet:
    a_init, _, c_init, _, _, _ = ucell2
    a, _, c, _, beta, _ = ucell
    a_ref, c_ref = RUC.ucell_manager.variables

    err_init = np.linalg.norm((abs(a_init-a)/a, abs(c_init-c)/c))*100
    err_ref = np.linalg.norm((abs(a_ref-a)/a, abs(c_ref-c)/c))*100

    print ("Percent error in unit cell before refinement: %2.7g %%" % err_init)
    print (" ''                ''      after refinement: %2.7g %%" % err_ref)
    assert err_ref < 1e-2*err_init
else:
    a_init, b_init, c_init, _, beta_init, _ = ucell2
    a, b, c, _, beta, _ = ucell
    a_ref, b_ref, c_ref, beta_ref = RUC.ucell_manager.variables
    beta_ref = beta_ref * 180 / np.pi

    err_init = np.linalg.norm((abs(a_init-a)/a, abs(b_init-b)/b, abs(c_init-c)/c, abs(beta_init-beta)/beta))*100
    err_ref = np.linalg.norm((abs(a_ref-a)/a, abs(b_ref-b)/b, abs(c_ref-c)/c, abs(beta_ref-beta)/beta))*100

    print ("\nPercent error in unit cell before refinement: %2.7g %%" % err_init)
    print (" ''                ''      after refinement: %2.7g %%\n" % err_ref)
    assert err_ref < 1e-2*err_init

print("OK!")
