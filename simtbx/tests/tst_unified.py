"""
Extend the previous tests in tst_gauss_argchk, tst_exafel_api, tst_kokkos_lib
New tests provide a context flag to exercise either cuda or kokkos_gpu
New tests implement a multipanel detector, not monolithic
Test 0) use monochromatic interface to generate a background-only pattern as reference
     1) exafel api interface, GPU background, monochromatic, use add_energy_channel_from_gpu_amplitudes
     2) exafel api interface, GPU background, monochromatic, use add_energy_multichannel_mask_allpanel (equate 1,2)
     3) exafel api interface, GPU background,  3-wavelength, use add_energy_channel_from_gpu_amplitudes
     4) exafel api interface, GPU background,  3-wavelength, use add_energy_multichannel_mask_allpanel (equate 3,4)
     5-6) verify the get_whitelist_raw_pixels API
Moved to tst_shoeboxes:
     7) implement shoebox masks and resulting mask
     8) 3-wavelength, use add_energy_channel_mask_allpanel, use all-pixel int mask (equate 3,8)
     9) 3-wavelength, use add_energy_channel_mask_allpanel, use all-pixel bool mask (equate 3,9)
Not implemented yet: multipanel detector, write to .h5 file
                     also implement with 3-color, not just mono, thus forcing weights
     10) 3-wavelength, use add_energy_multichannel_mask_allpanel, use all-pixel int mask (equate 3,10)
     11) 3-wavelength, use add_energy_multichannel_mask_allpanel, use whitelist
"""
from __future__ import absolute_import, division, print_function
import numpy as np
import dxtbx
from scitbx.array_family import flex
from dxtbx_model_ext import flex_Beam
from scitbx.matrix import sqr
from simtbx.nanoBragg import nanoBragg, shapetype
from simtbx.nanoBragg.tst_gauss_argchk import water, basic_crystal, basic_beam, basic_detector, amplitudes
from simtbx import get_exascale
from dxtbx.model.beam import BeamFactory
import math
from scitbx.math import five_number_summary

def parse_input():
  from iotbx.phil import parse
  master_phil="""
    context = kokkos_gpu *cuda
      .type = choice
      .optional = False
      .help = backend for parallel execution
  """
  phil_scope = parse(master_phil)
  # The script usage
  import libtbx.load_env # implicit import
  from dials.util.options import ArgumentParser
  # Create the parser
  parser = ArgumentParser(
        usage="\n libtbx.python tst_unified_api context=[kokkos_gpu|cuda]",
        phil=phil_scope,
        epilog="test monolithic detector, three-energy beam, cuda vs. kokkos")
  # Parse the command line. quick_parse is required for MPI compatibility
  params, options = parser.parse_args(show_diff_phil=True,quick_parse=True)
  return params,options

class several_wavelength_case:
  def __init__(self, BEAM, DETECTOR, CRYSTAL, SF_model, weights, special=0):
    SIM = nanoBragg(DETECTOR, BEAM, panel_id=0)
    nchannels = len(weights)
    assert nchannels in [1,3]
    if nchannels==1:
      print("\nassume monochromatic")
      self.wavlen = flex.double([BEAM.get_wavelength()])
    else:
      print("\nassume three energy channels")
      self.wavlen = flex.double([BEAM.get_wavelength()-0.02, BEAM.get_wavelength(), BEAM.get_wavelength()+0.02])
    self.frac = weights
    self.flux = self.frac*SIM.flux

    self.sfall_channels = {}
    for x in range(len(self.wavlen)):
      self.sfall_channels[x] = SF_model.get_amplitudes(at_angstrom = self.wavlen[x])
    self.DETECTOR = DETECTOR
    self.BEAM = BEAM
    self.CRYSTAL = CRYSTAL
    self.domains_per_crystal = 5.E10 # put Bragg spots on larger scale relative to background
    self.special = special # flag one-time special test cases

  def set_pythony_beams(self,SIM): # for the multiwavelength case with use of multiple nanoBragg sources
    pythony_beams = flex_Beam()
    for x in range(len(self.flux)):
      beam_descr = {'direction': (0.0, 0.0, 1.0),
             'divergence': 0.0,
             'flux': 1e11,
             'polarization_fraction': 1.,
             'polarization_normal': (0.0, 1.0, 0.0),
             'sigma_divergence': 0.0,
             'transmission': 1.0,
             'wavelength': self.wavlen[x]/1.e10} # not sure why this has to be in meters
      pythony_beams.append(BeamFactory.from_dict(beam_descr))
    SIM.xray_beams = pythony_beams

  def reset_pythony_beams(self,SIM): # for mono-wavelength addition of background
    pythony_beams = flex_Beam()
    beam_descr = {'direction': (0.0, 0.0, 1.0),
             'divergence': 0.0,
             'flux': 1e12,
             'polarization_fraction': 1.,
             'polarization_normal': (0.0, 1.0, 0.0),
             'sigma_divergence': 0.0,
             'transmission': 1.0,
             'wavelength': BEAM.get_wavelength()/1.e10} # not sure why this has to be in meters
    pythony_beams.append(BeamFactory.from_dict(beam_descr))
    SIM.xray_beams = pythony_beams

  def modularized_exafel_api_for_GPU(self, params, argchk=False, gpu_background=True, sources=False):
    gpu_channels_type = get_exascale("gpu_energy_channels",params.context)
    gpu_channels_singleton = gpu_channels_type (deviceId = 0)

    SIM = nanoBragg(self.DETECTOR, self.BEAM, panel_id=0)
    SIM.adc_offset_adu=0
    SIM.device_Id = 0

    assert gpu_channels_singleton.get_deviceID()==SIM.device_Id # API test
    assert gpu_channels_singleton.get_nchannels() == 0 # uninitialized
    for x in range(len(self.flux)):
          gpu_channels_singleton.structure_factors_to_GPU_direct(
           x, self.sfall_channels[x].indices(), self.sfall_channels[x].data())
    assert gpu_channels_singleton.get_nchannels() == len(self.flux) #API test
    SIM.Ncells_abc = (20,20,20)
    SIM.Amatrix = sqr(self.CRYSTAL.get_A()).transpose()
    SIM.oversample = 2
    SIM.xtal_shape = shapetype.Gauss
    SIM.interpolate = 0
    # allocate GPU arrays
    gpu_simulation = get_exascale("exascale_api",params.context)(nanoBragg = SIM)
    gpu_simulation.allocate()
    gpu_detector = get_exascale("gpu_detector",params.context)(
                 deviceId=SIM.device_Id, detector=self.DETECTOR, beam=self.BEAM)
    gpu_detector.each_image_allocate()
    # gpu_detector.show_summary()

    # two completely independent interfaces to loop over energies
    if sources is False: # original exascale API, explicit energy loop in Python
      for x in range(len(self.flux)):
        SIM.wavelength_A = self.wavlen[x]
        print("USE_EXASCALE_API+++++++++++++ Wavelength %d=%.6f, Flux %.6e, Fluence %.6e"%(
            x, SIM.wavelength_A, SIM.flux, SIM.fluence))
        gpu_simulation.add_energy_channel_from_gpu_amplitudes(
            x, gpu_channels_singleton, gpu_detector,
            weight = self.frac[x])
    else: # loop in C++; precludes using sources for divergence; uses sources from nanoBragg table
      self.set_pythony_beams(SIM)
      NN = 0 # compute the number of pixels
      for panel in self.DETECTOR:
          sz = panel.get_image_size()
          NN += sz[0]*sz[1]
      gpu_simulation.add_energy_multichannel_mask_allpanel(
            ichannels = flex.int(range(len(self.flux))),
            gpu_amplitudes = gpu_channels_singleton,
            gpu_detector = gpu_detector,
            pixel_active_list_ints = flex.size_t(range(NN)),
            weights = self.frac/len(self.frac)
      )
    per_image_scale_factor = self.domains_per_crystal # 1.0
    gpu_detector.scale_in_place(per_image_scale_factor) # apply scale directly on GPU
    if sources is False:
      SIM.wavelength_A = self.BEAM.get_wavelength() # return to canonical energy for subsequent background
    else:
      self.reset_pythony_beams(SIM)
    SIM.Fbg_vs_stol = water
    SIM.amorphous_sample_thick_mm = 0.02
    SIM.amorphous_density_gcm3 = 1
    SIM.amorphous_molecular_weight_Da = 18
    SIM.flux=1e12
    SIM.beamsize_mm=0.003 # square (not user specified)
    SIM.exposure_s=1.0 # multiplies flux x exposure
    gpu_simulation.add_background(gpu_detector)

    gpu_detector.write_raw_pixels(SIM)  # updates SIM.raw_pixels from GPU
    if self.special == 1: #special test breaks encapsulation
      self.special_test_case_1(SIM,gpu_detector)
      gpu_detector.each_image_free()
      return SIM
    self.data_array = gpu_detector.get_raw_pixels()
    assert self.data_array.focus()[0] == len(self.DETECTOR) # number of panels
    if len(self.DETECTOR) == 1: # if one panel, we can assert data are the same both ways
      single_size = self.data_array.focus()[1:3]
      self.data_array.reshape(flex.grid((single_size[0],single_size[1])))
      assert np.allclose(SIM.raw_pixels, self.data_array)
    # deallocate GPU arrays afterward
    gpu_detector.each_image_free()
    return SIM

  def special_test_case_1(self,SIM,gpu_detector):
    reference_raw_pixels = SIM.raw_pixels
    image_size = len(SIM.raw_pixels)
    even_pixels = flex.size_t(range(0,image_size,2))
    odd_pixels = flex.size_t(range(1,image_size,2))
    assert image_size == len(even_pixels) + len(odd_pixels)
    even_values = gpu_detector.get_whitelist_raw_pixels(even_pixels)
    odd_values = gpu_detector.get_whitelist_raw_pixels(odd_pixels)
    assert (len(even_values) == len(even_pixels)) and (len(odd_values) == len(odd_pixels))

    working_raw_pixels = flex.double(image_size) # blank array
    working_raw_pixels.set_selected(even_pixels, even_values)
    working_raw_pixels.set_selected(odd_pixels, odd_values)
    working_raw_pixels.reshape(flex.grid(SIM.raw_pixels.focus()))
    SIM.raw_pixels = working_raw_pixels
    SIM.to_cbf("test_unified_%s_005.cbf"%(params.context), intfile_scale=scale) # straight image; scale seems to come from global scope
    assert np.allclose(SIM.raw_pixels, reference_raw_pixels) # reassembly equates with original image

    workin2_raw_pixels = flex.double(image_size) # blank array
    workin2_raw_pixels.set_selected(even_pixels, odd_values) #assemble the image from two half-selections
    workin2_raw_pixels.set_selected(odd_pixels, even_values) #but in the wrong order
    workin2_raw_pixels.reshape(flex.grid(SIM.raw_pixels.focus()))
    SIM.raw_pixels = workin2_raw_pixels
    SIM.to_cbf("test_unified_%s_006.cbf"%(params.context), intfile_scale=scale) # mixed up checkerboard
    assert not np.allclose(SIM.raw_pixels, reference_raw_pixels) # wrong-order distorts image

if __name__=="__main__":
  params,options = parse_input()
  # make the dxtbx objects
  BEAM = basic_beam()
  DETECTOR = basic_detector()
  CRYSTAL = basic_crystal()
  SF_model = amplitudes(CRYSTAL)
  # Famp = SF_model.Famp # simple uniform amplitudes
  SF_model.random_structure(CRYSTAL)
  SF_model.ersatz_correct_to_P1()

  # Initialize based on GPU context
  gpu_instance_type = get_exascale("gpu_instance", params.context)
  gpu_instance = gpu_instance_type(deviceId = 0)

  print("\n# Use case 0 (%s). Background only"%params.context)
  SWC = several_wavelength_case(BEAM, DETECTOR, CRYSTAL, SF_model, weights=flex.double([0.]))
  SIM0 = SWC.modularized_exafel_api_for_GPU(params=params, argchk=False, gpu_background=True)

  print("\n# Use case 1 (%s). Monochromatic source"%params.context)
  SWC = several_wavelength_case(BEAM, DETECTOR, CRYSTAL, SF_model, weights=flex.double([1.]))
  SIM1 = SWC.modularized_exafel_api_for_GPU(params=params, argchk=False, gpu_background=True)
  scale = SIM1.get_intfile_scale() # case 1 is the overall scale factor to write all result files
  SIM1.to_cbf("test_unified_%s_001.cbf"%(params.context), intfile_scale=scale)
  SIM0.to_cbf("test_unified_%s_000.cbf"%(params.context), intfile_scale=scale)

  background_pixels = scale * SIM0.raw_pixels.as_1d()
  mean_background_pixel = flex.mean(background_pixels)
  print("mean",mean_background_pixel,"on",len(background_pixels))
  print("five number summary",five_number_summary(background_pixels))

  bragg_pixels_whole_image = (scale * (SIM1.raw_pixels - SIM0.raw_pixels).as_1d())
  bragg_pixels_only = bragg_pixels_whole_image.select(bragg_pixels_whole_image > 1.)
  print("mean",flex.mean(bragg_pixels_only),"on",len(bragg_pixels_only))
  print("five number summary",five_number_summary(bragg_pixels_only))
  # require the max Bragg peak to be 100x greater than the mean background
  assert flex.max(bragg_pixels_only) > 100. * mean_background_pixel

  multiple = 8.
  SWC_mult = several_wavelength_case(BEAM, DETECTOR, CRYSTAL, SF_model, weights=flex.double([multiple]))
  SIM1_mult = SWC_mult.modularized_exafel_api_for_GPU(params=params, argchk=False, gpu_background=True)
  # Bragg spot intensity should scale exactly with weight, using add_energy_channel_from_gpu_amplitudes
  assert np.allclose( (SIM1_mult.raw_pixels - SIM0.raw_pixels), multiple*(SIM1.raw_pixels - SIM0.raw_pixels) )

  if "kokkos" in params.context:
    print("\n# Use case 2 (%s). Monochromatic source"%params.context)#"sources=True": use multichannel API
    SIM2 = SWC.modularized_exafel_api_for_GPU(params=params, argchk=False, gpu_background=True, sources=True)
    SIM2.to_cbf("test_unified_%s_002.cbf"%(params.context), intfile_scale=scale)
  # Bragg spot intensity should scale exactly with weight, using multichannel API
    assert np.allclose(SIM1.raw_pixels, SIM2.raw_pixels) # equate monochromatic, two Bragg interfaces

  # Use explicitly 3-color code and equate with explicitly 1-color:
  print ("\n Equate background")
  weights = flex.double([(0./6.), (0./6.), (0./6.)]) # background only shot
  SWC = several_wavelength_case(BEAM, DETECTOR, CRYSTAL, SF_model, weights=weights)
  SIM_ctrl_py = SWC.modularized_exafel_api_for_GPU(params=params, argchk=False, gpu_background=True)
  assert np.allclose(SIM_ctrl_py.raw_pixels, SIM0.raw_pixels) # equate background
  if "kokkos" in params.context:
    SIM_ctrl_nb = SWC.modularized_exafel_api_for_GPU(params=params, argchk=False, gpu_background=True, sources=True)
    assert np.allclose(SIM_ctrl_nb.raw_pixels, SIM0.raw_pixels) # equate background

  print ("\n Equate 1-color")
  weights = flex.double([(0./6.), (6./6.), (0./6.)]) # implicitly one color
  SWC = several_wavelength_case(BEAM, DETECTOR, CRYSTAL, SF_model, weights=weights)
  SIM_ctrl_py = SWC.modularized_exafel_api_for_GPU(params=params, argchk=False, gpu_background=True)
  assert np.allclose(SIM_ctrl_py.raw_pixels, SIM1.raw_pixels) # equate one color
  if "kokkos" in params.context:
    SIM_ctrl_nb = SWC.modularized_exafel_api_for_GPU(params=params, argchk=False, gpu_background=True, sources=True)
    assert np.allclose(SIM_ctrl_nb.raw_pixels, SIM1.raw_pixels) # equate one color

  print ("\n Equate 2-color")
  weights = flex.double([(3./6.), (0./6.), (3./6.)]) # two color simulation with 3 explicit channels
  SWC = several_wavelength_case(BEAM, DETECTOR, CRYSTAL, SF_model, weights=weights)
  SIM_ctrl_py = SWC.modularized_exafel_api_for_GPU(params=params, argchk=False, gpu_background=True)
  SIM_ctrl_py.to_cbf("test_unified_%s_ctrl2_002.cbf"%(params.context), intfile_scale=scale)
  if "kokkos" in params.context:
    SIM_ctrl_nb = SWC.modularized_exafel_api_for_GPU(params=params, argchk=False, gpu_background=True, sources=True)
    assert np.allclose(SIM_ctrl_nb.raw_pixels, SIM_ctrl_py.raw_pixels) # equate two color with 2 different APIs

  # Controls
  # flux in 1:3:2 ratio
  # flux in 2:3:1 ratio
  print("\n# Use case 3 (%s). 3-Color"%params.context)
  weights = flex.double([(1./6.), (3./6.), (2./6.)]) # three color simulation slightly biased toward lower energy
  SWC = several_wavelength_case(BEAM, DETECTOR, CRYSTAL, SF_model, weights=weights)
  SIM3 = SWC.modularized_exafel_api_for_GPU(params=params, argchk=False, gpu_background=True)
  SIM3.to_cbf("test_unified_%s_003.cbf"%(params.context), intfile_scale=scale)
  if "kokkos" in params.context:
    print("\n# Use case 4 (%s). 3-Color"%params.context)
    SIM4 = SWC.modularized_exafel_api_for_GPU(params=params, argchk=False, gpu_background=True, sources=True)
    SIM4.to_cbf("test_unified_%s_004.cbf"%(params.context), intfile_scale=scale)
    assert np.allclose(SIM3.raw_pixels, SIM4.raw_pixels) # equate three color with 2 different APIs

  print("\n# Use case 3a (%s). 3-Color with bias toward high energy"%params.context)
  weights = flex.double([(2./6.), (3./6.), (1./6.)]) # three color simulation slightly biased toward higher energy
  SWC = several_wavelength_case(BEAM, DETECTOR, CRYSTAL, SF_model, weights=weights)
  SIM3a = SWC.modularized_exafel_api_for_GPU(params=params, argchk=False, gpu_background=True)
  SIM3a.to_cbf("test_unified_%s_103.cbf"%(params.context), intfile_scale=scale)

  # Now evaluate Bragg spots' radius of gyration (lower energy should be higher radius)
  low_energy_px = (scale*(SIM3.raw_pixels-SIM0.raw_pixels)).as_numpy_array()
  hi_energy_px = (scale*(SIM3a.raw_pixels-SIM0.raw_pixels)).as_numpy_array()
  xpts = np.array(range(low_energy_px.shape[0]))
  ypts = np.array(range(low_energy_px.shape[1]))
  X2D, Y2D = np.meshgrid(xpts,ypts)
  # center of mass = sum(mass * position) / sum(mass)
  center_of_mass_xy = (low_energy_px.shape[0]/2., low_energy_px.shape[1]/2.)
  # radius of gyration sq = sum(mass * (r_sq from com)) / sum(mass)
  lo_rg_sq = ( np.sum(low_energy_px * ((X2D-center_of_mass_xy[0])**2 + (Y2D-center_of_mass_xy[1])**2))/np.sum(low_energy_px)  )
  hi_rg_sq = ( np.sum(hi_energy_px * ((X2D-center_of_mass_xy[0])**2 + (Y2D-center_of_mass_xy[1])**2))/np.sum(hi_energy_px)  )
  print("low energy radius of gyration",math.sqrt(lo_rg_sq))
  print("hi energy radius of gyration",math.sqrt(hi_rg_sq))
  assert lo_rg_sq > hi_rg_sq

  print ("\n Use cases 5-6. Test the whitelist concept")
  # perform with 1-color image, same image as SIM1
  # recalculate the several_wavelength_case instance:
  reference = SIM1.raw_pixels
  SWC = several_wavelength_case(BEAM, DETECTOR, CRYSTAL, SF_model, weights=flex.double([1.]), special=1)
  SIM_working = SWC.modularized_exafel_api_for_GPU(params=params, argchk=False, gpu_background=True)

print("OK")
