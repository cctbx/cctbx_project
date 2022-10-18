"""
Extend the previous tests
Implemented unit tests for the whitelist mechanism:
     7) implement shoebox masks and resulting mask
     8) 3-wavelength, use add_energy_channel_mask_allpanel, use all-pixel int mask (equate 3,8)
     9) 3-wavelength, use add_energy_channel_mask_allpanel, use all-pixel bool mask (equate 3,9)
Unit test strategy:
Do a reference simulation with monochromatic beam.  Use dials to find strong spots.
Then construct two masks:  the positive mask just on the strong spot shoeboxes,
                       and the negative mask on all other pixels, including weak tails, weak spots, and background.
Then do two Bragg spot whitelist calls: one on tne positve mask, one on the negative mask.
A test image consisting of the sum of background + positive mask values + negative mask values should equal the reference image.
Furthermore, the result should be the same for either type of whitelist: bool array or integer addresses.
The GPU classes required modification for this to work, specifically the deepest array had to
be reinitialized after each Bragg calculation.
"""
from __future__ import absolute_import, division, print_function
import numpy as np
from dials.array_family import flex
from scitbx.matrix import sqr
from simtbx.tests.tst_unified import several_wavelength_case
from simtbx.nanoBragg.tst_gauss_argchk import water, basic_crystal, basic_beam, basic_detector, amplitudes
from simtbx import get_exascale
from simtbx.nanoBragg import nanoBragg, shapetype

def _dials_phil_str():
    return """
  include scope dials.algorithms.spot_finding.factory.phil_scope
"""
def spots_params():
  from iotbx.phil import parse
  dials_phil_str = _dials_phil_str()
  master_phil="""
    context = kokkos_gpu *cuda
      .type = choice
      .optional = False
      .help = backend for parallel execution
  output {
    shoeboxes = True
      .type = bool
      .help = Save the raw pixel values inside the reflection shoeboxes during spotfinding.
  }
  """
  phil_scope = parse(master_phil + dials_phil_str, process_includes=True)
  # The script usage
  import libtbx.load_env # implicit import
  from dials.util.options import ArgumentParser
  # Create the parser
  parser = ArgumentParser(usage="",phil=phil_scope, epilog="")
  # Parse the command line. quick_parse is required for MPI compatibility
  params, options = parser.parse_args(show_diff_phil=True,quick_parse=True)
  return params, options

class whitelist_case(several_wavelength_case):
  def extensions_for_shoebox_usage(self, params, argchk=False, gpu_background=True, sources=False):
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
    assert sources is False # original exascale API, explicit energy loop in Python

    for x in range(len(self.flux)):
        SIM.wavelength_A = self.wavlen[x]
        print("USE_WHITELIST(bool)_API+++++++++++++ Wavelength %d=%.6f, Flux %.6e, Fluence %.6e"%(
            x, SIM.wavelength_A, SIM.flux, SIM.fluence))
        gpu_simulation.add_energy_channel_mask_allpanel(
            x, gpu_channels_singleton, gpu_detector, positive_mask.as_1d())
            #weight = self.frac[x])
    per_image_scale_factor = self.domains_per_crystal # 1.0
    gpu_detector.scale_in_place(per_image_scale_factor) # apply scale directly on GPU
    positive_mask_bool_pixels = gpu_detector.get_raw_pixels()
    gpu_detector.scale_in_place(0) # reset
    """ explanation:  There are two arrays on device
    m_floatimage is the deepest one, that adds up bragg spots for each simulation call
    m_accumulate_floatimage is incremented by m_floatimage in the add_array(), only at the end of each simulation call
    the scale_in_place(0) only resets the shallow array (m_accumulate_floatimage)
    the deeper array (m_floatimage) could only be zeroed previously by the deallocate/allocate cycle
    instead, propose changing the add_array() kernel so that it does two jobs:
      1) add the deep array (m_floatimage) to the shallow (which is the current behavior)
      2) zero out the deep array (new behavior, should fix problem)
    """
    for x in range(len(self.flux)):
        SIM.wavelength_A = self.wavlen[x]
        print("USE_WHITELIST(bool)_API+++++++++++++ Wavelength %d=%.6f, Flux %.6e, Fluence %.6e"%(
            x, SIM.wavelength_A, SIM.flux, SIM.fluence))
        gpu_simulation.add_energy_channel_mask_allpanel(
            x, gpu_channels_singleton, gpu_detector, negative_mask.as_1d())
    gpu_detector.scale_in_place(per_image_scale_factor) # apply scale directly on GPU
    negative_mask_bool_pixels = gpu_detector.get_raw_pixels()
    gpu_detector.scale_in_place(0) # reset

    if True: # insert tests here, for the integer-address whitelist interface
      for x in range(len(self.flux)):
        SIM.wavelength_A = self.wavlen[x]
        print("USE_WHITELIST(pixel_address)_API+++++++++++++ Wavelength %d=%.6f, Flux %.6e, Fluence %.6e"%(
            x, SIM.wavelength_A, SIM.flux, SIM.fluence))
        gpu_simulation.add_energy_channel_mask_allpanel(
            x, gpu_channels_singleton, gpu_detector, positive_mask.iselection())
        per_image_scale_factor = self.domains_per_crystal # 1.0
        gpu_detector.scale_in_place(per_image_scale_factor) # apply scale directly on GPU
        positive_mask_int_pixels = gpu_detector.get_raw_pixels()
      gpu_detector.scale_in_place(0) # reset
      for x in range(len(self.flux)):
        SIM.wavelength_A = self.wavlen[x]
        print("USE_WHITELIST(pixel_address)_API+++++++++++++ Wavelength %d=%.6f, Flux %.6e, Fluence %.6e"%(
            x, SIM.wavelength_A, SIM.flux, SIM.fluence))
        gpu_simulation.add_energy_channel_mask_allpanel(
            x, gpu_channels_singleton, gpu_detector, negative_mask.iselection())
      gpu_detector.scale_in_place(per_image_scale_factor) # apply scale directly on GPU
      negative_mask_int_pixels = gpu_detector.get_raw_pixels()
      gpu_detector.scale_in_place(0) # reset
      assert np.allclose(positive_mask_int_pixels, positive_mask_bool_pixels) # either API gives same answer (int vs. bool)
      assert np.allclose(negative_mask_int_pixels, negative_mask_bool_pixels)

    SIM.wavelength_A = self.BEAM.get_wavelength() # return to canonical energy for subsequent background
    SIM.Fbg_vs_stol = water
    SIM.amorphous_sample_thick_mm = 0.02
    SIM.amorphous_density_gcm3 = 1
    SIM.amorphous_molecular_weight_Da = 18
    SIM.flux=1e12
    SIM.beamsize_mm=0.003 # square (not user specified)
    SIM.exposure_s=1.0 # multiplies flux x exposure
    gpu_simulation.add_background(gpu_detector)
    """Problem statement (before the GPU code was modified)
with pos first
pos -> back + pos
neg -> back + neg + pos
both-> back + neg + 2*pos

with neg first
pos -> back + neg + pos
neg -> back + neg
both-> back + 2*neg + pos

"""
    gpu_detector.write_raw_pixels(SIM)  # updates SIM.raw_pixels from GPU
    gpu_detector.each_image_free()
    SIM.raw_pixels += positive_mask_bool_pixels
    SIM.raw_pixels += negative_mask_bool_pixels
    return SIM

if __name__=="__main__":
  params,options = spots_params()
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

  print("\n# Use case 1 (%s). Monochromatic source"%params.context)
  SWC = several_wavelength_case(BEAM, DETECTOR, CRYSTAL, SF_model, weights=flex.double([1.]))
  SIM1 = SWC.modularized_exafel_api_for_GPU(params=params, argchk=False, gpu_background=True)
  scale = SIM1.get_intfile_scale() # case 1 is the overall scale factor to write all result files
  SIM1.to_cbf("test_shoebox_%s_001.cbf"%(params.context), intfile_scale=scale)
  SIM1.raw_pixels *= scale # temporarily set to intfile_scale for spotfinder
  from dxtbx.model.experiment_list import (
    Experiment,
    ExperimentList,
    ExperimentListFactory,
  )
  experiments = ExperimentListFactory.from_imageset_and_crystal(SIM1.imageset, CRYSTAL)

  # Find the strong spots
  observed = flex.reflection_table.from_observations(experiments, params, is_stills=True)
  observed.as_file("test_shoebox_%s_001.refl"%(params.context))
  SIM1.raw_pixels /= scale # reset from intfile_scale after spotfinder
  # Find the positive mask
  image_grid = flex.grid(SIM1.raw_pixels.focus())
  positive_mask = flex.bool(image_grid, False)

  imask = SIM1.imageset.get_mask(0)
  for x in range(len(observed)):
    fast1, fast2, slow1, slow2, z1, z2 = observed["bbox"][x]
    for slow in range(slow1,slow2):
      for fast in range(fast1,fast2):
        positive_mask[(slow,fast)]=True
  negative_mask=~positive_mask
  import pickle
  with open("test_shoebox_%s_001.mask"%(params.context),"wb") as M:
    pickle.dump([positive_mask], M)
  with open("test_shoebox_%s_002.mask"%(params.context),"wb") as M:
    pickle.dump([negative_mask], M)

  print("\n# Use case 2 (%s). Pixel masks."%params.context)
  print(positive_mask.count(True), negative_mask.count(True), len(positive_mask))
  WC = whitelist_case(BEAM, DETECTOR, CRYSTAL, SF_model, weights=flex.double([1.]))
  SIM3 = WC.extensions_for_shoebox_usage(params=params, argchk=False, gpu_background=True)
  SIM3.to_cbf("test_shoebox_%s_003.cbf"%(params.context), intfile_scale=scale)

  assert np.allclose(SIM1.raw_pixels, SIM3.raw_pixels) # reassembly equates with original image

print("OK")
