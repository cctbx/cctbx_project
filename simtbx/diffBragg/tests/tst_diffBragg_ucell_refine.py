from __future__ import division

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
parser.add_argument("--refine", action='store_true')
parser.add_argument("--new", action='store_true')
parser.add_argument("--crystalsystem", default='tetragonal',
                    choices=["monoclinic", "tetragonal"])
parser.add_argument("--curvatures", action='store_true')
args = parser.parse_args()
is_tet = args.crystalsystem == "tetragonal"
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.stats import pearsonr
import pylab as plt
from scipy.stats import linregress
from dxtbx.model import Experiment
from simtbx.nanoBragg import make_imageset
from simtbx.diffBragg.phil import phil_scope
from scitbx.matrix import sqr
from scitbx.matrix import rec

from simtbx.nanoBragg.nanoBragg_crystal import NBcrystal
from cctbx import uctbx
from simtbx.nanoBragg.sim_data import SimData
from simtbx.diffBragg import utils
from simtbx.diffBragg.refiners import crystal_systems
from dxtbx.model.crystal import Crystal

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
nbcryst = NBcrystal()
nbcryst.dxtbx_crystal = C
nbcryst.n_mos_domains = 1
nbcryst.thick_mm = 0.01
nbcryst.Ncells_abc = (7, 7, 7)

# STEP5: make an instance of diffBRagg, use the simData wrapper
SIM = SimData()
# overwrite the default detector with a smaller pixels one
img_sh = 700,700
SIM.detector = SimData.simple_detector(300, 0.1, img_sh)
SIM.crystal = nbcryst
SIM.instantiate_diffBragg(oversample=0, verbose=0, auto_set_spotscale=True)
# D is an instance of diffBragg with reasonable parameters
# and our dxtbx crystal created above
D = SIM.D
SIM.D.compute_curvatures = args.curvatures
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


# STEP7:
# compute the scattering and its derivative
print("Adding diffBragg spots")
#D.printout_pixel_fastslow =150, 351
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
shifts = [1e-4*(2*i) for i in range(1, 12, 2)]

if not args.refine:
    for i_param in range(n_ucell_params):
        analy_deriv = derivs[i_param]
        diffs = []
        diffs2 = []
        error = []
        error2 = []
        h_vals = []
        for i_shift, percent_shift in enumerate(shifts):

            if is_tet:
                var = [ucell[0], ucell[2]]
            else:
                var = [ucell[0], ucell[1], ucell[2], ucell[4]*np.pi/180]

            #name = UcellMan.variable_names[i_param]
            #if "alpha" in name or "beta" in name or "gamma" in name:
            #    param_shift = param_shift * 1e-4  # rescale the radian parameters

            param_shift = var[i_param] * percent_shift

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

            error.append(ave_error)
            h_vals.append( param_shift)
            print ("\tAverage error=%f; parameter shift h=%f" % (ave_error, abs(param_shift)))
            if args.curvatures:
                ave_error2 = np.abs(finite_second_deriv[bragg] - second_derivs[i_param][bragg]).mean()
                print("\tsecond derivative Average error=%f; parameter shift squared h^2=%f"
                      % (ave_error2, abs(param_shift)**2))
                r2 = pearsonr(second_derivs[i_param][bragg].ravel(), finite_second_deriv[bragg].ravel())[0]

                error2.append(ave_error2)
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

        l = linregress(h_vals, error)

        print ("finite diff l.rvalue=%10.7g" % l.rvalue)
        assert l.rvalue > .99
        assert l.slope > 0
        assert l.pvalue < 1e-6

        if args.curvatures:
            l2 = linregress(np.array(h_vals)**2, error2)
            print ("finite 2nd diff l.rvalue=%10.7g" % l2.rvalue)
            assert l2.rvalue > .99
            assert l2.slope > 0
            assert l2.pvalue < 1e-6

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
        assert np.all(np.diff(list(zip(trend[:-1], trend[1:])), axis=1) <= 0)
        if args.curvatures:
            assert (diffs2[0] > .99), "%f" % diffs2[0]
    print("OK!")

if is_tet:
    ucell2 = (55.5, 55.5, 76.1, 90, 90, 90)
else:  # args.crystalsystem == "monoclinic"
    ucell2 = (70.5, 60.1, 50.1, 90.0, 110.1, 90.0)

a2_real, b2_real, c2_real = sqr(uctbx.unit_cell(ucell2).orthogonalization_matrix()).transpose().as_list_of_lists()
C2 = Crystal(a2_real, b2_real, c2_real, symbol)

print ("Ensure the U matrices are the same for the ground truth and perturbation")
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
SIM.instantiate_diffBragg(oversample=0, auto_set_spotscale=True)
SIM.D.progress_meter = False
SIM.water_path_mm = 0.005
SIM.air_path_mm = 0.1
SIM.add_air = True
SIM.add_Water = True
SIM.include_noise = True
SIM.D.add_diffBragg_spots()
spots = SIM.D.raw_pixels.as_numpy_array()
SIM._add_background()
bg_img = SIM.D.raw_pixels.as_numpy_array()

SIM.D.readout_noise_adu = 3
SIM._add_noise()
# this is the ground truth image:
img = SIM.D.raw_pixels.as_numpy_array()
SIM.D.raw_pixels *= 0

E = Experiment()
E.detector = SIM.detector
E.beam = SIM.D.beam
E.crystal = C2  # intentionally set the wrong xtal model
E.imageset = make_imageset([img], E.beam, E.detector)

refls = utils.refls_from_sims([spots], E.detector, E.beam, thresh=20)

P = phil_scope.extract()
P.roi.shoebox_size = 20
P.roi.reject_edge_reflections = False
P.refiner.refine_Bmatrix = [1]
P.refiner.sensitivity.unitcell = [1, 1, 1, 0.1, 0.1, 0.1]
P.refiner.max_calls = [1000]
P.refiner.tradeps = 1e-10
# NOTE RUC.gtol = .9
# NOTE RUC.trad_conv = True  #False
# NOTE RUC.drop_conv_max_eps = 1e-9
P.refiner.curvatures = args.curvatures
P.refiner.use_curvatures_threshold = 0
P.refiner.poissononly = False
P.refiner.verbose = True
P.refiner.big_dump = False
P.refiner.sigma_r = SIM.D.readout_noise_adu
P.refiner.adu_per_photon = SIM.D.quantum_gain
P.simulator.crystal.has_isotropic_ncells = True
P.simulator.crystal.ncells_abc = 12, 12, 12
P.simulator.init_scale = SIM.D.spot_scale
P.simulator.beam.size_mm = SIM.beam.size_mm

# assert RUC.all_ang_off[0] < 0.005
from simtbx.diffBragg import refine_launcher
RUC = refine_launcher.local_refiner_from_parameters(refls, E, P, miller_data=SIM.crystal.miller_array)

if is_tet:
    a_init, _, c_init, _, _, _ = ucell2
    a, _, c, _, beta, _ = ucell
    a_ref, c_ref = RUC.UCELL_MAN[0].variables

    err_init = np.linalg.norm((abs(a_init-a)/a, abs(c_init-c)/c))*100
    err_ref = np.linalg.norm((abs(a_ref-a)/a, abs(c_ref-c)/c))*100

    print("Percent error in unit cell before refinement: %2.7g %%" % err_init)
    print(" ''                ''      after refinement: %2.7g %%" % err_ref)
    assert err_ref < 1e-2*err_init
else:
    a_init, b_init, c_init, _, beta_init, _ = ucell2
    a, b, c, _, beta, _ = ucell
    a_ref, b_ref, c_ref, beta_ref = RUC.UCELL_MAN[0].variables
    beta_ref = beta_ref * 180 / np.pi

    err_init = np.linalg.norm((abs(a_init-a)/a, abs(b_init-b)/b, abs(c_init-c)/c, abs(beta_init-beta)/beta))*100
    err_ref = np.linalg.norm((abs(a_ref-a)/a, abs(b_ref-b)/b, abs(c_ref-c)/c, abs(beta_ref-beta)/beta))*100

    print("\nPercent error in unit cell before refinement: %2.7g %%" % err_init)
    print(" ''                ''      after refinement: %2.7g %%\n" % err_ref)
    assert err_ref < 1e-2*err_init

print("OK!")
