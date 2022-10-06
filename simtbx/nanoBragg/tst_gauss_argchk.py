"""
Test the GAUSS_ARGCHK facility.  Basic idea:
The GAUSS shapetype uses a call to exp() measuring RLP distance to Ewald sphere.
Most pixels on a typical pattern are far from Bragg spots, so exp() evaluates to 0.
On GPU (but not CPU) we can save lots of execution time by pretesting the argument,
    for exp(-arg), evaluate to zero if arg >= 35
Provide a backdoor to the Mullen-Holton kernel, by defining shapetype=GAUSS_ARGCHK

This test exercises
1) the standard C++ result, i.e., the simple monochromatic case
2) the GPU implementation of the simple monochromatic case, if CUDA is enabled

The test is derived from tst_nanoBragg_cbf_write.py
Makes dxtbx models for detector, beam , crystal
Verifies pixel intensities are reproduced
"""
from __future__ import absolute_import, division, print_function
import numpy as np
from scipy import constants
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal

from cctbx import sgtbx, miller
from cctbx.crystal import symmetry

from dxtbx.model.beam import BeamFactory
from dxtbx.model.crystal import CrystalFactory
from dxtbx.model.detector import DetectorFactory
from scitbx.matrix import sqr, col
from simtbx.nanoBragg import nanoBragg, shapetype

# rough approximation to water: interpolation points for sin(theta/lambda) vs structure factor
water = flex.vec2_double([(0,2.57),(0.0365,2.58),(0.07,2.8),(0.12,5),(0.162,8),(0.18,7.32),(0.2,6.75),(0.216,6.75),(0.236,6.5),(0.28,4.5),(0.3,4.3),(0.345,4.36),(0.436,3.77),(0.5,3.17)])

def basic_crystal():
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
  return CrystalFactory.from_dict(cryst_descr)

def basic_beam():
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
  return BeamFactory.from_dict(beam_descr)

def basic_detector():
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
  return DetectorFactory.from_dict(det_descr)

class amplitudes:
  def __init__(self, CRYSTAL):
    # make a dummy HKL table with constant HKL intensity
    # this is just to make spots
    DEFAULT_F = 1e2
    symbol = CRYSTAL.get_space_group().info().type().lookup_symbol()  # this is just P43212
    assert symbol == "P 43 21 2" # test case, start with P43212, make P1 for nanoBragg
    sgi = sgtbx.space_group_info(symbol)
    symm = symmetry(unit_cell=CRYSTAL.get_unit_cell(), space_group_info=sgi)
    miller_set = symm.build_miller_set(anomalous_flag=True, d_min=1.6, d_max=999)
    Famp = flex.double(np.ones(len(miller_set.indices())) * DEFAULT_F)
    self.Famp = miller.array(miller_set=miller_set, data=Famp).set_observation_type_xray_amplitude()

  def random_structure(self,crystal):
    """We're going to do some very approximate stuff here.  Given a unit
     cell & SG, will put typical atomic contents in the unit cell & get
     structure factors.
    """
    import random
    random.seed(0)
    from scitbx.array_family import flex
    flex.set_random_seed(0)
    from cctbx.development import random_structure

    uc_volume = crystal.get_unit_cell().volume()
    asu_volume = uc_volume / crystal.get_space_group().order_z()
    target_number_scatterers = int(asu_volume)//128 # Very approximate rule of thumb for proteins with ~50% solvent content
    element_unit = ['O']*19 + ['N']*18 + ['C']*62 + ['S']*1 + ['Fe']*1
    element_pallet = element_unit * (1 + ( target_number_scatterers//len(element_unit) ))
    assert len(element_pallet) >= target_number_scatterers
    # Ersatz hard limit to prevent excessive execution time of xray_structure() below.
    elements = element_pallet[:min(1000, target_number_scatterers)]

    xs = random_structure.xray_structure(
      space_group_info = crystal.get_space_group().info(), unit_cell = crystal.get_unit_cell(),
      elements=elements, min_distance=1.2)
    self.xs = xs

  def ersatz_correct_to_P1(self):
    primitive_xray_structure = self.xs.primitive_setting()
    P1_primitive_xray_structure = primitive_xray_structure.expand_to_p1()
    self.xs = P1_primitive_xray_structure

  def get_amplitudes(self, at_angstrom):
    # Since we are getting amplitudes for nanoBragg, let us assure they are in P1
    symbol = self.xs.space_group().info().type().lookup_symbol()
    assert symbol=="P 1", "Must be in P1 to accept amplitudes for ExaFEL GPU interface"
    # take a detour to insist on calculating anomalous contribution of every atom
    scatterers = self.xs.scatterers()
    for sc in scatterers:
      from cctbx.eltbx import henke
      expected_henke = henke.table(sc.element_symbol()).at_angstrom(at_angstrom)
      sc.fp = expected_henke.fp()
      sc.fdp = expected_henke.fdp()

    import mmtbx.command_line.fmodel
    phil2 = mmtbx.command_line.fmodel.fmodel_from_xray_structure_master_params
    params2 = phil2.extract()
    params2.high_resolution = 1.6
    params2.fmodel.k_sol = 0.35
    params2.fmodel.b_sol = 46.
    params2.structure_factors_accuracy.algorithm = "fft"
    params2.output.type = "real"
    import mmtbx
    f_model = mmtbx.utils.fmodel_from_xray_structure(
      xray_structure = self.xs,
      f_obs          = None,
      add_sigmas     = True,
      params         = params2).f_model
    #f_model.show_summary()
    return f_model

CPU_GPU_Lookup = dict(add_nanoBragg_spots="CPU", add_nanoBragg_spots_cuda="GPU")

def simple_monochromatic_case(bragg_engine, BEAM, DETECTOR, CRYSTAL, SF_model, argchk=False):
  Famp = SF_model.get_amplitudes(at_angstrom=BEAM.get_wavelength())

  # do the simulation
  SIM = nanoBragg(DETECTOR, BEAM, panel_id=0)
  SIM.Ncells_abc = (20,20,20)
  SIM.Fhkl = Famp
  SIM.Amatrix = sqr(CRYSTAL.get_A()).transpose()
  SIM.oversample = 2
  if argchk:
    print("\nmonochromatic case,",CPU_GPU_Lookup[bragg_engine.__name__],"argchk")
    SIM.xtal_shape = shapetype.Gauss_argchk
  else:
    print("\nmonochromatic case,",CPU_GPU_Lookup[bragg_engine.__name__],"no argchk")
    SIM.xtal_shape = shapetype.Gauss
  bragg_engine(SIM) # appropriate add_nanoBragg_spots, either CPU or GPU
  domains_per_crystal = 5.E10 # put Bragg spots on larger scale relative to background
  SIM.raw_pixels*=domains_per_crystal
  ref_max_bragg = flex.max(SIM.raw_pixels) # get the maximum pixel value for a Bragg spot

  SIM.Fbg_vs_stol = water
  SIM.amorphous_sample_thick_mm = 0.02
  SIM.amorphous_density_gcm3 = 1
  SIM.amorphous_molecular_weight_Da = 18
  SIM.flux=1e12
  SIM.beamsize_mm=0.003 # square (not user specified)
  SIM.exposure_s=1.0 # multiplies flux x exposure
  SIM.progress_meter=False
  SIM.add_background()
  ref_mean_with_background = flex.mean(SIM.raw_pixels)
  print ("Ratio",ref_max_bragg/ref_mean_with_background)
  assert ref_max_bragg > 10. * ref_mean_with_background # data must be sensible, Bragg >> solvent
  return SIM

if __name__=="__main__":
  # make the dxtbx objects
  BEAM = basic_beam()
  DETECTOR = basic_detector()
  CRYSTAL = basic_crystal()
  SF_model = amplitudes(CRYSTAL)
  # Famp = SF_model.Famp # simple uniform amplitudes
  SF_model.random_structure(CRYSTAL)
  SF_model.ersatz_correct_to_P1()
  import sys
  runmode = sys.argv[1]
  assert runmode in ["CPU","GPU"]

  # Use case 1.  Simple monochromatic X-rays
  bragg_engine = nanoBragg.add_nanoBragg_spots
  SIM = simple_monochromatic_case(bragg_engine, BEAM, DETECTOR, CRYSTAL, SF_model, argchk=False)
  SIM2 = simple_monochromatic_case(bragg_engine, BEAM, DETECTOR, CRYSTAL, SF_model, argchk=True)
  output_scale = SIM.get_intfile_scale(intfile_scale=0)
  #print("The output scale is",output_scale)
  SIM.adc_offset_adu=0
  SIM.to_smv_format(fileout="test_full_001.img", intfile_scale=output_scale)
  assert approx_equal(SIM.raw_pixels, SIM2.raw_pixels)
  SIM.to_cbf("test_full_001.cbf", intfile_scale=output_scale)

  if runmode=="GPU":
    bragg_engine = nanoBragg.add_nanoBragg_spots_cuda
    SIM3 = simple_monochromatic_case(bragg_engine, BEAM, DETECTOR, CRYSTAL, SF_model, argchk=False)
    SIM3.to_cbf("test_full_003.cbf", intfile_scale=output_scale)
    SIM4 = simple_monochromatic_case(bragg_engine, BEAM, DETECTOR, CRYSTAL, SF_model, argchk=True)
    assert approx_equal(SIM.raw_pixels, SIM3.raw_pixels)
    assert approx_equal(SIM.raw_pixels, SIM4.raw_pixels)

  print("OK")
