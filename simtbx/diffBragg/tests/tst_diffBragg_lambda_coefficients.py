
"""
This test checks the lambda coefficients property and derivatives
"""

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action='store_true')
parser.add_argument("--idx", type=int, help="coefficient index (0 or 1)", default=0, choices=[0,1])
args = parser.parse_args()

import numpy as np
import pylab as plt
from scipy.stats import linregress
from scipy.spatial.transform import Rotation
from simtbx.diffBragg import sim_data
from scitbx.matrix import sqr, rec
from cctbx import uctbx
from dxtbx.model import Crystal

ucell = (70, 60, 50, 90.0, 110, 90.0)
symbol = "C121"

a_real, b_real, c_real = sqr(uctbx.unit_cell(ucell).orthogonalization_matrix()).transpose().as_list_of_lists()
C = Crystal(a_real, b_real, c_real, symbol)

# random raotation
rotation = Rotation.random(num=1, random_state=101)[0]
Q = rec(rotation.as_quat(), n=(4, 1))
rot_ang, rot_axis = Q.unit_quaternion_as_axis_and_angle()
C.rotate_around_origin(rot_axis, rot_ang)

S = sim_data.SimData()
S.crystal.dxtbx_crystal = C
spectrum = S.beam.spectrum
wave, flux = spectrum[0]
Nwave = 5
waves = np.linspace(wave-wave*0.002, wave+wave*0.002, Nwave)
fluxes = np.ones(Nwave) * flux / Nwave

lambda0_GT = waves[0]
lambda1_GT = waves[1] - waves[0]
assert np.allclose(waves, np.arange(Nwave) * lambda1_GT + lambda0_GT)

S.beam.spectrum = list(zip(waves, fluxes))
S.detector = sim_data.SimData.simple_detector(180, 0.1, (1024, 1024))
S.instantiate_diffBragg(verbose=0, oversample=0)
S.D.lambda_coefficients = lambda0_GT, lambda1_GT
S.D.spot_scale = 100000
S.D.Ncells_abc = 12

S.D.refine(12)
S.D.initialize_managers()
S.D.region_of_interest = ((0, 1023), (0, 1023))

S.D.add_diffBragg_spots()
img = S.D.raw_pixels.as_numpy_array()
derivs = S.D.get_lambda_derivative_pixels()
deriv0 = derivs[0].as_numpy_array()
deriv1 = derivs[1].as_numpy_array()

S.D.raw_pixels *= 0
S.D.use_lambda_coefficients = False
S.D.add_diffBragg_spots()
test_img = S.D.raw_pixels.as_numpy_array()
assert np.allclose(img,test_img)
S.D.use_lambda_coefficients = True
S.D.raw_pixels *= 0
print ("OK")

bragg = img > 1e-1  # select bragg scattering regions

all_error = []
all_error2 = []
shifts = []
shifts2 = []

from scipy import constants
ENERGY_CONV = 1e10*constants.c*constants.h / constants.electron_volt

energy_shifts = 0.1, .3, .5, 1, 3, 5, 10   # in electron volt
b_percs = 0.001, 0.002, 0.004, 0.008,  0.016, 0.032, 0.064
reference_energy = ENERGY_CONV / wave
for i_shift, en_shift in enumerate(energy_shifts):

    wave_shifted = ENERGY_CONV / (reference_energy + en_shift)
    wave_shift = wave - wave_shifted
    delta_a = wave_shift

    delta_b = lambda1_GT*b_percs[i_shift]

    if args.idx == 0:
        new_waves = np.arange(Nwave)*(lambda1_GT) + lambda0_GT + delta_a
        shift = delta_a
    else:
        new_waves = np.arange(Nwave) * (lambda1_GT + delta_b) + lambda0_GT
        shift = delta_b

    en = np.mean(ENERGY_CONV/new_waves)

    if args.idx == 0:
        S.D.lambda_coefficients = lambda0_GT + shift, lambda1_GT
        shifts.append(shift)
    else:
        S.D.lambda_coefficients = lambda0_GT, lambda1_GT + shift
        shifts.append(shift)

    S.D.raw_pixels *= 0
    S.D.region_of_interest = ((0, 1023), (0, 1023))
    S.D.add_diffBragg_spots()
    img2 = S.D.raw_pixels.as_numpy_array()

    fdiff = (img2 - img) / shift

    if args.idx == 0:
        error = np.abs(fdiff[bragg] - deriv0[bragg]).mean()
    else:
        error = np.abs(fdiff[bragg] - deriv1[bragg]).mean()

    all_error.append(error)

    print ("error=%f, step=%f, energy=%f" % (error, delta_a, en))
    #if args.plot:
    #    plt.subplot(121)
    #    plt.imshow(fdiff)
    #    plt.title("finite diff")
    #    plt.subplot(122)
    #    plt.imshow(deriv)
    #    plt.title("analytical")
    #    plt.draw()
    #    plt.suptitle("Shift %d / %d"
    #                 % (i_shift + 1, len(perc)))
    #    plt.pause(0.8)

if args.plot:
    #plt.close()
    plt.plot(shifts, all_error, 'o')
    plt.show()
    if args.curvatures:
        plt.plot(shifts2, all_error2, 'o')
        plt.show()

l = linregress(shifts, all_error)
assert l.rvalue > .9999  # this is definitely a line!
assert l.slope > 0
assert l.pvalue < 1e-6
#if args.curvatures:
#    l = linregress(shifts2, all_error2)
#    assert l.rvalue > .9999  # this is definitely a line!
#    assert l.slope > 0
#    assert l.pvalue < 1e-6

print("OK!")
