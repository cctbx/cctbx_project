from __future__ import division

from dxtbx.model.crystal import Crystal
from cctbx import uctbx
from scitbx.matrix import sqr, rec, col
import numpy as np
from scipy.spatial.transform import Rotation
from scitbx.matrix import sqr

from simtbx.nanoBragg.nanoBragg_crystal import NBcrystal
from simtbx.nanoBragg.sim_data import SimData
from simtbx.diffBragg.utils import fcalc_from_pdb
from simtbx.nanoBragg.utils import ENERGY_CONV
from scipy.signal import windows


ucell = (79, 79, 38, 90,90,90)
symbol = "P43212"

miller_array_GT = fcalc_from_pdb(resolution=2, wavelength=1, algorithm='fft', ucell=ucell, symbol=symbol)
Ncells_gt = 12, 12, 12

np.random.seed(3142019)
# generate a random rotation
rotation = Rotation.random(num=1, random_state=100)[0]
Q = rec(rotation.as_quat(), n=(4, 1))
rot_ang, rot_axis = Q.unit_quaternion_as_axis_and_angle()

a_real, b_real, c_real = sqr(uctbx.unit_cell(ucell).orthogonalization_matrix()).transpose().as_list_of_lists()
x = col((-1, 0, 0))
y = col((0, -1, 0))
z = col((0, 0, -1))
rx, ry, rz = np.random.uniform(-180, 180, 3)
RX = x.axis_and_angle_as_r3_rotation_matrix(rx, deg=True)
RY = y.axis_and_angle_as_r3_rotation_matrix(ry, deg=True)
RZ = z.axis_and_angle_as_r3_rotation_matrix(rz, deg=True)
M = RX * RY * RZ
a_real = M * col(a_real)
b_real = M * col(b_real)
c_real = M * col(c_real)
C = Crystal(a_real, b_real, c_real, symbol)
C.rotate_around_origin(rot_axis, rot_ang)

# Setup the simulation and create a realistic image
# with background and noise
# <><><><><><><><><><><><><><><><><><><><><><><><><>
nbcryst = NBcrystal(init_defaults=True)
nbcryst.dxtbx_crystal = C   # simulate ground truth
nbcryst.thick_mm = 0.1
nbcryst.Ncells_abc = Ncells_gt  # ground truth Ncells

nbcryst.miller_array = miller_array_GT
print("Ground truth ncells = %f" % (nbcryst.Ncells_abc[0]))

# ground truth detector
DET_gt = SimData.simple_detector(150, 0.177, (600, 600))

# initialize the simulator
SIM = SimData(use_default_crystal=True)
spec = SIM.beam.spectrum
total_flux = spec[0][1]
wave = spec[0][0]

en = ENERGY_CONV / wave
delta_en = 1.5
ens_gt = np.arange(en - 5, en + 6, delta_en)
waves_gt = ENERGY_CONV / ens_gt
num_energies = len(ens_gt)
fluxes_gt = np.ones(num_energies) * total_flux / num_energies
fluxes_gt = fluxes_gt*windows.hann(num_energies)
fluxes_gt /= fluxes_gt.sum()
fluxes_gt *= total_flux

spectrum_GT = list(zip(waves_gt, fluxes_gt))
gt_lambda0 = waves_gt[0]
gt_lambda1 = waves_gt[1] - waves_gt[0]
spec_idx = np.arange(num_energies)
assert np.allclose(waves_gt, gt_lambda0 + spec_idx*gt_lambda1)

np.random.seed(3142019)
lam0 = np.random.normal(gt_lambda0, gt_lambda0 * 0.002)
lam1 = np.random.normal(gt_lambda1, abs(gt_lambda1) * 0.002)
waves_perturbed = lam0 + spec_idx * lam1
print("ENERGY TRUTH=%.4f" % (ENERGY_CONV / gt_lambda0))
print("ENERGY PERTURBED=%.4f" % (ENERGY_CONV / lam0))

SIM.beam.spectrum = list(zip(waves_gt, fluxes_gt))
SIM.detector = DET_gt
SIM.crystal = nbcryst
SIM.instantiate_diffBragg(oversample=0, auto_set_spotscale=True)
SIM.D.nopolar = False
SIM.D.default_F = 0
SIM.D.progress_meter = False
SIM.water_path_mm = 0.15
SIM.air_path_mm = 0.1
SIM.add_air = True
SIM.add_water = True
SIM.include_noise = True
SIM.D.add_diffBragg_spots()
SPOTS = SIM.D.raw_pixels.as_numpy_array()
SIM.D.readout_noise_adu = 1
SIM._add_background()
SIM._add_noise()

# This is the ground truth image:
img = SIM.D.raw_pixels.as_numpy_array()
SIM.D.raw_pixels *= 0

from simtbx.diffBragg import utils
from dxtbx.model import Experiment
from simtbx.nanoBragg import make_imageset
from cctbx_project.simtbx.diffBragg.phil import phil_scope
from simtbx.diffBragg import refine_launcher

E = Experiment()
E.detector = DET_gt
E.beam = SIM.beam.nanoBragg_constructor_beam
E.crystal = C
E.imageset = make_imageset([img], E.beam, E.detector)

refls = utils.refls_from_sims([SPOTS], E.detector, E.beam, thresh=20)
spectra_file = "_testing_spectra_file.lam"
np.savetxt(spectra_file, np.vstack((waves_perturbed, fluxes_gt)).T, delimiter=',', header="wavelengths, weights")

P = phil_scope.extract()
P.roi.shoebox_size = 20
P.roi.reject_edge_reflections = False
P.refiner.refine_spectra = [1]
P.refiner.sensitivity.spectra_coefficients = 1, 1
P.refiner.init.spectra_coefficients = 0, 1
P.refiner.ranges.spectra0 = -0.01, 0.01
P.refiner.ranges.spectra1 = 0.95, 1.05
P.simulator.spectrum.filename = spectra_file
P.simulator.beam.size_mm = SIM.beam.size_mm
P.simulator.total_flux = total_flux
P.refiner.max_calls = [20]
P.refiner.tradeps = 1e-7
P.refiner.curvatures = False
P.refiner.use_curvatures_threshold = 0
P.refiner.poissononly = False
P.refiner.verbose = True
P.refiner.big_dump = False
P.refiner.sigma_r = SIM.D.readout_noise_adu
P.refiner.adu_per_photon = SIM.D.quantum_gain
P.simulator.crystal.has_isotropic_ncells = True
P.simulator.init_scale = SIM.D.spot_scale
P.simulator.crystal.ncells_abc = Ncells_gt
P.refiner.ncells_mask = "111"

RUC = refine_launcher.local_refiner_from_parameters(refls, E, P, miller_data=SIM.crystal.miller_array)
coef = RUC._get_spectra_coefficients()
waves_refined = coef[0] + coef[1] * waves_perturbed
fluxsum = sum(fluxes_gt)
en_ref_com = ENERGY_CONV / (sum(fluxes_gt * waves_refined) / fluxsum)
en_com = ENERGY_CONV / (sum(fluxes_gt * waves_gt) / fluxsum)
en_init_com = ENERGY_CONV / (sum(fluxes_gt * waves_perturbed) / fluxsum)

print("Before refinement: COM energy=%f" % en_init_com)
print("AFTER refinement: COM energy=%f" % en_ref_com)
print("Ground truth COM energy = %f" % en_com)
assert abs(en_ref_com - en_com) < 1

print("OK")
