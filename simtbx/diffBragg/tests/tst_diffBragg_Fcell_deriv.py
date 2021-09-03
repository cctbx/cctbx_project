from __future__ import division

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
parser.add_argument("--curvatures", action='store_true')
args = parser.parse_args()

if args.plot:
    import pylab as plt
    plt.figure()
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.stats import linregress

from scitbx.array_family import flex
from scitbx.matrix import sqr
from scitbx.matrix import rec

from simtbx.nanoBragg.nanoBragg_crystal import NBcrystal
from cctbx import uctbx
from simtbx.nanoBragg.sim_data import SimData
from dxtbx.model.crystal import Crystal

ucell = (70, 60, 50, 90.0, 110, 90.0)
symbol = "C121"

a_real, b_real, c_real = sqr(uctbx.unit_cell(ucell).orthogonalization_matrix()).transpose().as_list_of_lists()
C = Crystal(a_real, b_real, c_real, symbol)

# random raotation
rotation = Rotation.random(num=1, random_state=101)[0]
Q = rec(rotation.as_quat(), n=(4, 1))
rot_ang, rot_axis = Q.unit_quaternion_as_axis_and_angle()
C.rotate_around_origin(rot_axis, rot_ang)

# make a nanoBragg crystal to pass to diffBragg
nbcryst = NBcrystal(init_defaults=True)
nbcryst.dxtbx_crystal = C
nbcryst.n_mos_domains = 1
nbcryst.thick_mm = 0.01
nbcryst.Ncells_abc = (7, 7, 7)

# make an instance of diffBRagg, use the simData wrapper
SIM = SimData(use_default_crystal=True)
# overwrite the default detector with a smaller pixels one
SIM.detector = SimData.simple_detector(300, 0.1, (700, 700))
SIM.crystal = nbcryst
Fcell = 1e6
SIM.instantiate_diffBragg(oversample=0, verbose=0, interpolate=0, default_F=Fcell)
# D is an instance of diffBragg with reasonable parameters
# and our dxtbx crystal created above
D = SIM.D

D.progress_meter = True

# initialize the derivative manager for Fcell
fcell = 11  # internal index of fcell manager within diffBragg
D.refine(fcell)
D.initialize_managers()

#roi = ((0, 699), (0, 699))
#rX = slice(roi[0][0], roi[0][1], 1)
#rY = slice(roi[1][0], roi[1][1], 1)
#D.region_of_interest = roi

# compute the scattering and its derivative
print("Adding diffBragg spots")
# Set all miller indices to have the same Fcell value
indices, data = D.Fhkl_tuple
data = flex.double(np.ones(len(indices))*Fcell)
D.F000 = Fcell
D.Fhkl_tuple = indices, data.deep_copy(), None
D.compute_curvatures = True

# NOTE optionally focus on a single Bragg peak
#x1,x2 = 335, 365
#y1,y2 = 155, 185
#D.region_of_interest = (x1, x2), (y1, y2)
D.add_diffBragg_spots()

print("Done!")
img = D.raw_pixels_roi.as_numpy_array()

# reset all pixel values
D.raw_pixels *= 0
D.raw_pixels_roi *= 0
first_deriv = D.get_derivative_pixels(fcell).as_numpy_array()
if args.curvatures:
    second_deriv = D.get_second_derivative_pixels(fcell).as_numpy_array()

# iterate over the parameters and do a finite difference test for each one
shifts = [1e-2*(i**2) for i in range(1, 20)]
#shifts = np.array([.5,1,2,4,8,16,32]) # 1e-3, 1e-2, 1e-1, 1, 10]
all_error = []
all_delta_F = []

for i_shift, percent_shift in enumerate(shifts):

    delta_F = Fcell * percent_shift*1e-2

    D.F000 = Fcell + delta_F
    D.Fhkl_tuple = indices, data.deep_copy()+delta_F, None
    D.default_F = Fcell+delta_F
    D.add_diffBragg_spots()

    img_forward = D.raw_pixels_roi.as_numpy_array()
    # reset for next computation
    D.raw_pixels_roi *= 0
    D.raw_pixels *= 0

    finite_deriv = (img_forward-img) / delta_F
    ave_error = np.abs(finite_deriv - first_deriv).mean()
    #ave_error = np.abs(finite_deriv[bragg] - first_deriv[bragg]).mean()
    all_error.append(ave_error)
    all_delta_F.append(delta_F)

    if args.curvatures:
        # estimate the second derivative
        D.F000 = Fcell - delta_F
        D.Fhkl_tuple = indices, data.deep_copy() - delta_F, None
        D.default_F = Fcell - delta_F
        D.add_diffBragg_spots()
        img_backward = D.raw_pixels_roi.as_numpy_array()

        # reset for next computation
        D.raw_pixels_roi *= 0
        D.raw_pixels *= 0

        finite_second_deriv = (img_forward - 2*img + img_backward) / delta_F / delta_F
        assert np.allclose(finite_second_deriv, second_deriv)
        # NOTE: second derivative is a constant, doesnt depend on F_cell, hence these are always super close

    print("error=%2.7g, step=%f (%d/%d)" % (ave_error, delta_F, i_shift + 1, len(shifts)))

if args.plot:
    plt.close()
    plt.plot(all_delta_F, all_error, 'o')
    plt.show()

# one expects errors to scale linearly with parameter shift for finite differences
l = linregress(all_delta_F, all_error)
assert l.rvalue > .9999  # this is definitely a line!
assert l.slope > 0
assert l.pvalue < 1e-6
print("Error versus parameter shift fits a line with slope=%2.7g and Correleation Coef=%2.7g" % (l.slope, l.rvalue))
print("OK!")
