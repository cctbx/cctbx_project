from __future__ import division
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
parser.add_argument("--crystalsystem", default='tetragonal',
                    choices=["monoclinic", "tetragonal"])
parser.add_argument("--curvatures", action='store_true')
parser.add_argument("--cuda", action="store_true")
args = parser.parse_args()

if args.cuda:
    import os
    os.environ["DIFFBRAGG_USE_CUDA"]="1"

import numpy as np
from scipy.spatial.transform import Rotation
from scipy.stats import pearsonr
import pylab as plt
from scipy.stats import linregress
from simtbx.diffBragg.hopper_utils import UCELL_ID_OFFSET
from scitbx.matrix import sqr
from scitbx.matrix import rec
from simtbx.nanoBragg.nanoBragg_crystal import NBcrystal
from cctbx import uctbx
from simtbx.nanoBragg.sim_data import SimData
from simtbx.diffBragg import utils
from dxtbx.model.crystal import Crystal


# STEP 1:
# make a crystal and orient it randomly
if args.crystalsystem=="tetragonal":
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

# STEP3:
# create the unit cell parameter manager
UcellMan = utils.manager_from_params(ucell)
n_ucell_params = len(UcellMan.variables)
assert np.allclose(UcellMan.B_recipspace, C.get_B())

# STEP4:
# make a nanoBragg crystal to pass to diffBragg
nbcryst = NBcrystal(init_defaults=True)
nbcryst.dxtbx_crystal = C
nbcryst.n_mos_domains = 1
nbcryst.thick_mm = 0.01
nbcryst.Ncells_abc = (7, 7, 7)

# STEP5: make an instance of diffBRagg, use the simData wrapper
SIM = SimData()
# overwrite the default detector to use smaller pixels
img_sh = 700,700
SIM.detector = SimData.simple_detector(300, 0.1, img_sh)
SIM.crystal = nbcryst
SIM.instantiate_diffBragg(oversample=0, verbose=0, auto_set_spotscale=True)
# D is an instance of diffBragg with sensible parameters
# and our dxtbx crystal created above
D = SIM.D
D.progress_meter = True

# STEP6:
# initialize the derivative managers for the unit cell parameters
for i_param in range(n_ucell_params):
    D.refine(UCELL_ID_OFFSET+i_param)
for i in range(n_ucell_params):
    D.set_ucell_derivative_matrix(UCELL_ID_OFFSET+i, UcellMan.derivative_matrices[i])
    if args.curvatures:
        D.set_ucell_second_derivative_matrix(UCELL_ID_OFFSET + i, UcellMan.second_derivative_matrices[i])
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
    analy_deriv = D.get_derivative_pixels(UCELL_ID_OFFSET+i_param).as_numpy_array()
    derivs.append(analy_deriv)
    if args.curvatures:
        second_derivs.append(D.get_second_derivative_pixels(UCELL_ID_OFFSET+i_param).as_numpy_array())

# STEP8
# iterate over the parameters and do a finite difference test for each one
# parameter shifts:
shifts = [1e-4*(2*i) for i in range(1, 12, 2)]

import copy
starting_var = copy.copy(UcellMan.variables)
for i_param in range(n_ucell_params):
    analy_deriv = derivs[i_param]
    cc_vals = []
    cc_vals2 = []
    error = []
    error2 = []
    h_vals = []
    for i_shift, percent_shift in enumerate(shifts):

        var = copy.copy(starting_var)

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
        cc_vals.append(r)

        error.append(ave_error)
        h_vals.append( param_shift)
        print ("\tAverage error=%f; parameter shift h=%f" % (ave_error, abs(param_shift)))
        if args.curvatures:
            ave_error2 = np.abs(finite_second_deriv[bragg] - second_derivs[i_param][bragg]).mean()
            print("\tsecond derivative Average error=%f; parameter shift squared h^2=%f"
                  % (ave_error2, abs(param_shift)**2))
            r2 = pearsonr(second_derivs[i_param][bragg].ravel(), finite_second_deriv[bragg].ravel())[0]

            error2.append(ave_error2)
            cc_vals2.append(r2)
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
        plt.plot(shifts, cc_vals, 'o')
        title = "Unit cell parameter %d / %d" % (i_param+1, n_ucell_params)
        plt.title(title + "\nPearson corr between finite deriv and analytical")
        plt.xlabel("unit cell shifts")
        plt.ylabel("Pearson corr")
        plt.show()
        if args.curvatures:
            plt.close()
            plt.plot(np.array(shifts)**2, cc_vals2, 'o')
            title = "Unit cell parameter %d / %d" % (i_param + 1, n_ucell_params)
            plt.title(title + "\nPearson corr between finite second deriv and analytical")
            plt.xlabel("unit cell shifts")
            plt.ylabel("Pearson corr")
            plt.show()

    # verify a high correlation for the smallest parameter shift
    print("Check high pearson R between analytical and finite diff")
    print("Pearson correlection at smallest parameter shift=%f" % cc_vals[0])
    assert(cc_vals[0] > .98), "%f" % cc_vals[0]
    # check monotonic decrease
    print("Fit polynomial and check monotonic decrease")
    trend = np.polyval(np.polyfit(shifts, cc_vals, 2), shifts)
    assert np.all(np.diff(list(zip(trend[:-1], trend[1:])), axis=1) <= 0)
    if args.curvatures:
        assert (cc_vals2[0] > .99), "%f" % cc_vals2[0]
print("OK!")
