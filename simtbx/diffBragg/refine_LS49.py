
import os
datdir = "/Users/dermen/crystal/modules/LS49/tests/ls49_big_data"
os.environ["LS49_BIG_DATA"] = datdir
from copy import deepcopy

import dxtbx
from scitbx.matrix import sqr
from scitbx.array_family import flex
from dxtbx.model import Crystal

from LS49.spectra.generate_spectra import spectra_simulation
from LS49.sim.util_fmodel import gen_fmodel
from LS49.sim.step5_pad import data
from LS49.sim.step4_pad import microcrystal

from simtbx.diffBragg import helpers
from simtbx.diffBragg.sim_data2 import SimData
from simtbx.diffBragg.nanoBragg_crystal import nanoBragg_crystal
from simtbx.diffBragg.nanoBragg_beam import nanoBragg_beam
from simtbx.diffBragg.load_ls49 import process_ls49_image

out = process_ls49_image()
spotroi = out["bboxes_x1x2y1y2"]
idx_crystal = out["dxcrystal"]
abc_init = out["tilt_abc"]
data_img = out["data_img"]
spectrum = out["spectrum"]

poly = False

loader = dxtbx.load(os.path.join(datdir, "ls49_0.npz"))
ls49_det = loader.get_detector()
ls49_beam = loader.get_beam()

# STEP 1: get a spectrum
SS = spectra_simulation()
I = SS.generate_recast_renormalized_images(20, energy=7120, total_flux=1e12)
out = next(I)
wavelens = out[0]
fluxes = out[1]
wavelen_A = out[2]

if not poly:
    wavelens = [wavelen_A]
    fluxes = [sum(fluxes)]

# STEP 2 , make the structure factors
local_data = data()
GF = gen_fmodel(resolution=1.7,
                pdb_text=local_data.get("pdb_lines"),
                algorithm="fft",
                wavelength=wavelen_A)
GF.set_k_sol(0.435)
GF.make_P1_primitive()
sfall_main = GF.get_amplitudes()
sg_symbol = sfall_main.space_group_info().type().lookup_symbol()
ucell_tuple = sfall_main.unit_cell().parameters()

# make the crystal parameters
crystal = microcrystal(Deff_A=4000, length_um=4, beam_diameter_um=1)
N = crystal.number_of_cells(sfall_main.unit_cell())
print("Ncells abc is %d|%d|%d" % (N,N,N))

# unique orientation
mt = flex.mersenne_twister(seed=0)
rotation = sqr(flex.mersenne_twister(seed=0).random_double_r3_rotation_matrix())

# make the nanoBragg crystal
C = nanoBragg_crystal()
C.Ncells_abc = (N, N, N)
C.mos_spread_deg = 0.05
C.n_mos_domains = 25

# grab the ground truth a,b,c
C.dxtbx_crystal = nanoBragg_crystal.dxtbx_crystal_from_ucell_and_symbol(ucell_tuple, "P1")  #indexing_result_C
C.missetting_matrix = rotation
a_gt,b_gt,c_gt = C.a_b_c_realspace_misset

# set the indexing result
C.dxtbx_crystal = idx_crystal
C.missetting_matrix = sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))

spec_idx = 0
GF.reset_wavelength(wavelens[spec_idx])
GF.reset_specific_at_wavelength(
    label_has="FE1", tables=local_data.get("Fe_oxidized_model"), newvalue=wavelens[spec_idx])
GF.reset_specific_at_wavelength(
    label_has="FE2", tables=local_data.get("Fe_reduced_model"), newvalue=wavelens[spec_idx])
print("USING scatterer-specific energy-dependent scattering factors")
sfall_channel = GF.get_amplitudes()
C.miller_array = sfall_channel

# make the beam
B = nanoBragg_beam()
B.spectrum = spectrum  # zip(wavelens, fluxes)
B.size_mm = 0.003
B.unit_s0 = ls49_beam.get_unit_s0()

S = SimData()
S.beam = B
S.crystal = C
S.detector = ls49_det
S.seed = 1
S.water_path_mm = 0.1
S.air_path_mm = 10

brargs = {
    'adc_offset': 10,
    'interpolate': 0,
    'default_F': 0,
    'verbose': 0,
    'oversample': 1}
S.instantiate_diffBragg(**brargs)
S.update_nanoBragg_instance("spot_scale", crystal.domains_per_crystal)
print("spot scale is %f" % crystal.domains_per_crystal)
S.D.show_params()

from simtbx.diffBragg.refiners import RefineRot
RR = RefineRot(spot_rois=spotroi, abc_init=abc_init,
               img=data_img, SimData_instance=S, plot_images=False)

# Step 2:
RXYZ = RR.get_rotation_misset()  # perturbation rotation matrix

a_init, b_init, c_init = C.a_b_c_realspace_misset
C_init = Crystal(a_init, b_init, c_init, "P1")

C_refined = deepcopy(C_init)
q = RXYZ.r3_rotation_matrix_as_unit_quaternion()
angle_RXYZ, axis_RXYZ = q.unit_quaternion_as_axis_and_angle(deg=True)
C_refined.rotate_around_origin(axis_RXYZ, angle_RXYZ)

angles = helpers.compare_with_ground_truth(
    a_gt, b_gt, c_gt, [C_init, C_refined])

print ("\nInitial indexing result misset: %f" % angles[0])
print ("Optimized misset: %f" % angles[1])

assert angles[1] < angles[0]
print("OK!")

#from simtbx.diffBragg import utils
#angles_rad = RR.x[-4:-1]
#C_refined = utils.refine_model_from_angles(C_init, angles=angles_rad)

#crystal.dxtbx_crystal_with_missetting()
##C.set_A(sqr(crystal.Amatrix_realspace).inverse())
#angles_rad = RR.x[-4:-1]
#from dials.algorithms.indexing import compare_orientation_matrices
#ang = compare_orientation_matrices.difference_rotation_matrix_axis_angle(C_init, C0)[2]
#ang2 = compare_orientation_matrices.difference_rotation_matrix_axis_angle(C_refined, C0)[2]
#
#from IPython import embed
#embed()

#assert ang > 0.5
#assert ang2 < 0.002
#print("OK!")


#S.add_water = True
#S.add_air = True
#S.include_noise = True
#
#S._add_nanoBragg_spots()
#reference = np.load("ls49_0_mono_bragg.npy")
#from IPython import embed
#embed()
#
#
#S._add_background()
#reference = np.load("ls49_0_mono_noiseless.npy")
#from IPython import embed
#embed()
#
#S._add_noise()
#img = S.D.raw_pixels.as_numpy_array()
#reference = np.load("ls49_0_mono.npy")
#embed()

