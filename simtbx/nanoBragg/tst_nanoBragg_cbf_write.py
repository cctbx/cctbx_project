"""
Makes dxtbx models for detector, beam , crystal

Uses models to instantiate nanoBragg and compute Bragg spots

Writes results to a full CBF using new nanoBragg method to_cbf (in nanoBragg/__init__.py)

Loads CBF with dxtbx, and uses the loaded detector and beam to recompute the Bragg spots

Verifies pixel intensities are reproduced
"""
from __future__ import absolute_import, division, print_function
import numpy as np
from scipy import constants

from cctbx import sgtbx, miller
from cctbx.crystal import symmetry
import dxtbx
from dxtbx.model.beam import BeamFactory
from dxtbx.model.crystal import CrystalFactory
from dxtbx.model.detector import DetectorFactory
from scitbx.array_family import flex
from scitbx.matrix import sqr, col
from simtbx.nanoBragg import nanoBragg, shapetype

print("Make a randomly oriented xtal")
# make a randomly oriented crystal..
np.random.seed(3142019)
# make random rotation about principle axes
x = col((-1, 0, 0))
y = col((0, -1, 0))
z = col((0, 0, -1))
rx, ry, rz = np.random.uniform(-180, 180, 3)
RX = x.axis_and_angle_as_r3_rotation_matrix(rx, deg=True)
RY = y.axis_and_angle_as_r3_rotation_matrix(ry, deg=True)
RZ = z.axis_and_angle_as_r3_rotation_matrix(rz, deg=True)
M = RX*RY*RZ
real_a = M*col((79, 0, 0))
real_b = M*col((0, 79, 0))
real_c = M*col((0, 0, 38))
# dxtbx crystal description
cryst_descr = {'__id__': 'crystal',
               'real_space_a': real_a.elems,
               'real_space_b': real_b.elems,
               'real_space_c': real_c.elems,
               'space_group_hall_symbol': ' P 4nw 2abw'}

print("Make a beam")
# make a beam
ENERGY = 9000
ENERGY_CONV = 1e10*constants.c*constants.h / constants.electron_volt
WAVELEN = ENERGY_CONV/ENERGY
# dxtbx beam model description
beam_descr = {'direction': (0.0, 0.0, 1.0),
             'divergence': 0.0,
             'flux': 1e11,
             'polarization_fraction': 1.,
             'polarization_normal': (0.0, 1.0, 0.0),
             'sigma_divergence': 0.0,
             'transmission': 1.0,
             'wavelength': WAVELEN}

# make a detector panel
# monolithic camera description
print("Make a dxtbx detector")
detdist = 100.
pixsize = 0.1
im_shape = 1536, 1536
det_descr = {'panels':
               [{'fast_axis': (1.0, 0.0, 0.0),
                 'slow_axis': (0.0, -1.0, 0.0),
                 'gain': 1.0,
                 'identifier': '',
                 'image_size': im_shape,
                 'mask': [],
                 'material': '',
                 'mu': 0.0,
                 'name': 'Panel',
                 'origin': (-im_shape[0]*pixsize/2., im_shape[1]*pixsize/2., -detdist),
                 'pedestal': 0.0,
                 'pixel_size': (pixsize, pixsize),
                 'px_mm_strategy': {'type': 'SimplePxMmStrategy'},
                 'raw_image_offset': (0, 0),
                 'thickness': 0.0,
                 'trusted_range': (-1e7, 1e7),
                 'type': ''}]}

# make the dxtbx objects
BEAM = BeamFactory.from_dict(beam_descr)
DETECTOR = DetectorFactory.from_dict(det_descr)
CRYSTAL = CrystalFactory.from_dict(cryst_descr)

# make a dummie HKL table with constant HKL intensity
# this is just to make spots
DEFAULT_F = 1e2
symbol = CRYSTAL.get_space_group().info().type().lookup_symbol()  # this is just P43212
sgi = sgtbx.space_group_info(symbol)
symm = symmetry(unit_cell=CRYSTAL.get_unit_cell(), space_group_info=sgi)
miller_set = symm.build_miller_set(anomalous_flag=True, d_min=1.6, d_max=999)
Famp = flex.double(np.ones(len(miller_set.indices())) * DEFAULT_F)
Famp = miller.array(miller_set=miller_set, data=Famp).set_observation_type_xray_amplitude()

Ncells_abc = 20, 20, 20
oversmaple = 2

# do the simulation
print("Do the initial simulation")
SIM = nanoBragg(DETECTOR, BEAM, panel_id=0)
SIM.Ncells_abc = Ncells_abc
SIM.Fhkl = Famp
SIM.Amatrix = sqr(CRYSTAL.get_A()).transpose()
SIM.oversample = oversmaple
SIM.xtal_shape = shapetype.Gauss
SIM.add_nanoBragg_spots()

# Take note that the CBF writer and SMV writer handle spot scale differently
output_scale=1.E9
SIM.adc_offset_adu=0
smv_filename = "test_full_100.img"
SIM.to_smv_format(fileout=smv_filename, intfile_scale=output_scale)

# write the simulation to disk using cbf writer
cbf_filename = "test_full_100.cbf"
print("write simulation to disk (%s)" % cbf_filename)
SIM.to_cbf(cbf_filename, intfile_scale=output_scale)

# load the CBF from disk
print("Open file %s using dxtbx" % cbf_filename)
loader = dxtbx.load(cbf_filename)

print("Perform second simulation using dxtbx loaded detector and beam")
test_cbf = nanoBragg(loader.get_detector(), loader.get_beam(), panel_id=0)
test_cbf.Ncells_abc = Ncells_abc
test_cbf.Fhkl = Famp
# FIXME: for some reason cbf_writer converts the beam polarization fraction to 0.999,
# which would break this test hence we must set it back to 1..
test_cbf.polarization = 1
test_cbf.Amatrix = sqr(CRYSTAL.get_A()).transpose()
test_cbf.oversample = oversmaple
test_cbf.xtal_shape = shapetype.Gauss
test_cbf.add_nanoBragg_spots()

# verify test_cbf and SIM produce the same Bragg spot image
print("Check the intensities haven't changed, and that  cbf writing preserved geometry")
assert np.allclose(SIM.raw_pixels.as_numpy_array(), test_cbf.raw_pixels.as_numpy_array())

# verify cbf (double) and smv (int) produce the same image to within an ADU
loader_smv = dxtbx.load(smv_filename)
assert np.allclose(loader.get_raw_data().as_numpy_array(), loader_smv.get_raw_data().as_numpy_array(), atol=1.1)

print("OK!")
