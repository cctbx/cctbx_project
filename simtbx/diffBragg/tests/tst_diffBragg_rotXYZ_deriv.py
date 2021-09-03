from __future__ import division

from argparse import ArgumentParser
parser = ArgumentParser("diffBragg tests")
parser.add_argument("--plot", action='store_true')
parser.add_argument("--cuda", action='store_true')
parser.add_argument("--curvatures", action='store_true')
parser.add_argument("--plotimages", action="store_true")
parser.add_argument("--rotidx", type=int, choices=[0,1,2], required=True)
parser.add_argument("--randomrotate", type=int, default=None, help="seed for random rotation of Umatrix")
args = parser.parse_args()
if args.plot:
    import pylab as plt
from simtbx.nanoBragg.sim_data import SimData

from simtbx.nanoBragg.nanoBragg_crystal import NBcrystal
import numpy as np
from scitbx.matrix import sqr
from cctbx import uctbx
from dxtbx.model import Crystal
from scipy import stats

rot_idx = args.rotidx

ucell = (70, 60, 50, 90.0, 110, 90.0)
symbol = "C121"

a_real, b_real, c_real = sqr(uctbx.unit_cell(ucell).orthogonalization_matrix()).transpose().as_list_of_lists()
C = Crystal(a_real, b_real, c_real, symbol)
# make a nanoBragg crystal to pass to diffBragg
nbcryst = NBcrystal()
nbcryst.dxtbx_crystal = C
nbcryst.n_mos_domains = 1
nbcryst.thick_mm = 0.01
nbcryst.Ncells_abc = (7, 7, 7)

# make an instance of diffBRagg, use the simData wrapper
SIM = SimData(use_default_crystal=True)
# overwrite the default detector with a smaller pixels one
SIM.detector = SimData.simple_detector(220, 0.1, (1000, 1000))
SIM.crystal = nbcryst

SIM.instantiate_diffBragg(oversample=0, verbose=0, interpolate=0, default_F=1e3,auto_set_spotscale=True)
# D is an instance of diffBragg with reasonable parameters
# and our dxtbx crystal created above
D = SIM.D
if args.curvatures:
    D.compute_curvatures = True

# STEP 1: simulate the un-perturbed image:
D.refine(rot_idx)
D.initialize_managers()
D.set_value(rot_idx, 0)
D.use_cuda = args.cuda
#D.printout_pixel_fastslow = 786, 567
D.add_diffBragg_spots()
img0 = D.raw_pixels_roi.as_numpy_array()
bragg = img0 > 1  #np.ones_like(img0).astype(bool)
if args.plotimages:
    plt.title("Scattering from crystal")
    plt.imshow(img0)
    plt.show()

deriv = D.get_derivative_pixels(rot_idx).as_numpy_array()
if args.curvatures:
    second_deriv = D.get_second_derivative_pixels(rot_idx).as_numpy_array()

error_vals = []
error_vals2 = []
delta_h = []
delta_h2 = []
theta_vals = [0.0005 + i*0.0005 for i in range(8)]   # 0.01, 0.02048, 0.04096, 0.08192, 0.16384, 0.32768, 0.65536, 1.31072, 2.62144

for theta_degrees in theta_vals:
    theta = theta_degrees * np.pi / 180
    D.set_value(rot_idx, theta)

    # simulate the scattering in the rotated crystal:
    D.raw_pixels_roi *= 0
    D.add_diffBragg_spots()
    img_plus = D.raw_pixels_roi.as_numpy_array()
    if args.curvatures:
        D.set_value(rot_idx, -theta)
        D.raw_pixels_roi *= 0
        D.add_diffBragg_spots()
        img_minus = D.raw_pixels_roi.as_numpy_array()
        finite_second_diff = (img_plus-2*img0 + img_minus) / theta / theta
        delta_h2.append(theta**2)

    # STEP3 : compute finite differenceL
    finite_diff = (img_plus-img0)/theta
    delta_h.append(theta)

    if args.plotimages:
        plt.subplot(121)
        plt.imshow(finite_diff)
        plt.title("finite diff.")
        plt.subplot(122)
        plt.imshow(deriv)
        plt.title("analytical")
        plt.suptitle("Theta = %f deg. " % theta_degrees)
        plt.draw()
        plt.pause(1)
        if args.curvatures:
            plt.clf()
            plt.subplot(121)
            plt.imshow(finite_second_diff)
            plt.title("finite 2nd diff.")
            plt.subplot(122)
            plt.imshow(second_deriv)
            plt.title("analytical")
            plt.suptitle("Theta = %f deg. " % theta_degrees)
            plt.draw()
            plt.pause(0.2)

    error_image = abs(deriv[bragg]-finite_diff[bragg])
    error = error_image.mean()
    error_vals.append(error)

    if args.curvatures:
        error_image2 = abs(second_deriv[bragg]-finite_second_diff[bragg])
        error2 = error_image2.mean()
        print("Theta = %.4f deg, 1st error = %2.7g, 2nd %2.7g" %
                (theta_degrees, error, error2))
        error_vals2.append(error2)
    else:
        print("Theta = %.4f deg, 1st error = %2.7g" %
              (theta_degrees, error))

if args.plot:
    plt.close()
    plt.plot(delta_h, error_vals, '.')
    ax = plt.gca()
    ax.set_xlabel("theta (rad.)")
    ax.set_ylabel(r"$\langle |\,$finite_diff - analytical$\,| \rangle$", fontsize=14)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.subplots_adjust(left=0.12)
    plt.show()
    if args.curvatures:
        plt.plot(delta_h2, error_vals2, '.')
        ax = plt.gca()
        ax.set_xlabel("theta^2 (rad.^2)")
        ax.set_ylabel(r"$\langle |\,$finite_second_diff - analytical$\,| \rangle$", fontsize=14)
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        plt.subplots_adjust(left=0.12)
        plt.show()

l = stats.linregress(delta_h, error_vals)
assert l.rvalue > .9999, "%2.7g" % l.rvalue
assert l.slope > 0, "%2.7g" % l.slope
assert l.pvalue < 1e-6, "%2.7g" % l.pvalue
if args.curvatures:
    l = stats.linregress(delta_h2, error_vals2)
    assert l.rvalue > .9999, "2nd deriv rvalue %2.7g" % l.rvalue
    assert l.slope > 0, "2nd deriv slope %2.7g" % l.slope
    assert l.pvalue < 1e-6, "2nd deriv pvalue %2.7g" % l.pvalue
if args.cuda:
    D.gpu_free()
print("OK!")
