"""
Extend the tests in tst_gauss_argchk.
Both test cases excercise a simple monolithic detector
tst_exafel_api introduces polychromatic beam (3 wavelengths in this simplified example)
Supply command-line context=[kokkos_gpu | cuda] to switch between both implementations
Test 1) standard C++ result CPU
     2) CPU, but use the stable sort version of CPU background (FUTURE PLAN)
     3) exafel api interface to GPU, fast evaluation on many energy channels, CPU background
     4) exafel api interface to GPU, fast evaluation on many energy channels, GPU background
     5) exafel api interface to GPU, with GAUSS_ARGCHK
"""
from __future__ import absolute_import, division, print_function
import numpy as np
import dxtbx
from scitbx.array_family import flex
from scitbx.matrix import sqr
from simtbx.nanoBragg import nanoBragg, shapetype
from simtbx.nanoBragg.tst_gauss_argchk import water, basic_crystal, basic_beam, basic_detector, amplitudes
from simtbx import get_exascale

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
        usage="\n libtbx.python tst_exafel_api context=[kokkos_gpu|cuda]",
        phil=phil_scope,
        epilog="test monolithic detector, three-energy beam, cuda vs. kokkos")
  # Parse the command line. quick_parse is required for MPI compatibility
  params, options = parser.parse_args(show_diff_phil=True,quick_parse=True)
  return params,options

class several_wavelength_case:
 def __init__(self, BEAM, DETECTOR, CRYSTAL, SF_model):
  SIM = nanoBragg(DETECTOR, BEAM, panel_id=0)
  print("\nassume three energy channels")
  self.wavlen = flex.double([BEAM.get_wavelength()-0.002, BEAM.get_wavelength(), BEAM.get_wavelength()+0.002])
  self.flux = flex.double([(1./6.)*SIM.flux, (3./6.)*SIM.flux, (2./6.)*SIM.flux])
  self.sfall_channels = {}
  for x in range(len(self.wavlen)):
    self.sfall_channels[x] = SF_model.get_amplitudes(at_angstrom = self.wavlen[x])
  self.DETECTOR = DETECTOR
  self.BEAM = BEAM
  self.CRYSTAL = CRYSTAL
  self.domains_per_crystal = 5.E10 # put Bragg spots on larger scale relative to background

 def several_wavelength_case_for_CPU(self):
  SIM = nanoBragg(self.DETECTOR, self.BEAM, panel_id=0)
  SIM.adc_offset_adu=0
  for x in range(len(self.wavlen)):
    SIM.flux = self.flux[x]
    SIM.wavelength_A = self.wavlen[x]
    print("CPUnanoBragg_API+++++++++++++ Wavelength %d=%.6f, Flux %.6e, Fluence %.6e"%(
            x, SIM.wavelength_A, SIM.flux, SIM.fluence))
    SIM.Fhkl = self.sfall_channels[x]
    SIM.Ncells_abc = (20,20,20)
    SIM.Amatrix = sqr(self.CRYSTAL.get_A()).transpose()
    SIM.oversample = 2
    SIM.xtal_shape = shapetype.Gauss
    SIM.interpolate = 0
    SIM.add_nanoBragg_spots()
  SIM.raw_pixels*=self.domains_per_crystal
  ref_max_bragg = flex.max(SIM.raw_pixels) # get the maximum pixel value for a Bragg spot

  SIM.wavelength_A = self.BEAM.get_wavelength()
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

 def modularized_exafel_api_for_GPU(self, params, argchk=False, gpu_background=True):
  gpu_channels_type = get_exascale("gpu_energy_channels",params.context)
  gpu_channels_singleton = gpu_channels_type (deviceId = 0)

  SIM = nanoBragg(self.DETECTOR, self.BEAM, panel_id=0)
  SIM.adc_offset_adu=0
  SIM.device_Id = 0

  assert gpu_channels_singleton.get_deviceID()==SIM.device_Id
  assert gpu_channels_singleton.get_nchannels() == 0 # uninitialized
  for x in range(len(self.flux)):
          gpu_channels_singleton.structure_factors_to_GPU_direct(
           x, self.sfall_channels[x].indices(), self.sfall_channels[x].data())
  assert gpu_channels_singleton.get_nchannels() == len(self.flux)
  SIM.Ncells_abc = (20,20,20)
  SIM.Amatrix = sqr(self.CRYSTAL.get_A()).transpose()
  SIM.oversample = 2
  if argchk:
    print("\npolychromatic GPU argchk")
    SIM.xtal_shape = shapetype.Gauss_argchk
  else:
    print("\npolychromatic GPU no argchk")
    SIM.xtal_shape = shapetype.Gauss
  SIM.interpolate = 0
  # allocate GPU arrays
  gpu_simulation = get_exascale("exascale_api",params.context)(nanoBragg = SIM)
  gpu_simulation.allocate()

  gpu_detector = get_exascale("gpu_detector",params.context)(
                 deviceId=SIM.device_Id, detector=self.DETECTOR, beam=self.BEAM)
  gpu_detector.each_image_allocate()

  # loop over energies
  for x in range(len(self.flux)):
      SIM.flux = self.flux[x]
      SIM.wavelength_A = self.wavlen[x]
      print("USE_EXASCALE_API+++++++++++++ Wavelength %d=%.6f, Flux %.6e, Fluence %.6e"%(
            x, SIM.wavelength_A, SIM.flux, SIM.fluence))
      gpu_simulation.add_energy_channel_from_gpu_amplitudes(
        x, gpu_channels_singleton, gpu_detector)
  per_image_scale_factor = self.domains_per_crystal # 1.0
  gpu_detector.scale_in_place(per_image_scale_factor) # apply scale directly on GPU
  SIM.wavelength_A = self.BEAM.get_wavelength() # return to canonical energy for subsequent background

  if gpu_background:
      SIM.Fbg_vs_stol = water
      SIM.amorphous_sample_thick_mm = 0.02
      SIM.amorphous_density_gcm3 = 1
      SIM.amorphous_molecular_weight_Da = 18
      SIM.flux=1e12
      SIM.beamsize_mm=0.003 # square (not user specified)
      SIM.exposure_s=1.0 # multiplies flux x exposure
      gpu_simulation.add_background(gpu_detector)

      # deallocate GPU arrays afterward
      gpu_detector.write_raw_pixels(SIM)  # updates SIM.raw_pixels from GPU
      gpu_detector.each_image_free()
  else:
      # deallocate GPU arrays up front
      gpu_detector.write_raw_pixels(SIM)  # updates SIM.raw_pixels from GPU
      gpu_detector.each_image_free()

      SIM.Fbg_vs_stol = water
      SIM.amorphous_sample_thick_mm = 0.02
      SIM.amorphous_density_gcm3 = 1
      SIM.amorphous_molecular_weight_Da = 18
      SIM.flux=1e12
      SIM.beamsize_mm=0.003 # square (not user specified)
      SIM.exposure_s=1.0 # multiplies flux x exposure
      SIM.progress_meter=False
      SIM.add_background()
  return SIM

def diffs(labelA, A, labelB, B):
  diff = A-B
  min = flex.min(diff); mean = flex.mean(diff); max = flex.max(diff)
  print("Pixel differences between %s and %s, minimum=%.4f mean=%.4f maximum=%.4f"%(
       labelA, labelB, min, mean, max))
  assert min > -1.0
  assert max < 1.0

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

  print("\n# Use case 2.  Three-wavelength polychromatic source")
  SWC = several_wavelength_case(BEAM, DETECTOR, CRYSTAL, SF_model)
  SIM = SWC.several_wavelength_case_for_CPU()
  SIM.to_smv_format(fileout="test_full_cpu_002.img") # scales by default
  scale = SIM.get_intfile_scale()
  print ("Scale",scale)
  SIM.to_cbf("test_full_cpu_002.cbf", intfile_scale=scale)
  # verify cbf (double) and smv (int) produce the same image to within an ADU
  loader_smv = dxtbx.load("test_full_cpu_002.img")
  loader_cbf = dxtbx.load("test_full_cpu_002.cbf")
  assert np.allclose(loader_cbf.get_raw_data().as_numpy_array(), loader_smv.get_raw_data().as_numpy_array(), atol=1.1)

  # Switch the remaining tests based on GPU context
  gpu_instance_type = get_exascale("gpu_instance", params.context)
  gpu_instance = gpu_instance_type(deviceId = 0)

  print("\n# Use case 3 (%s): modularized api argchk=False, gpu_background=False"%params.context)
  SIM3 = SWC.modularized_exafel_api_for_GPU(params=params, argchk=False, gpu_background=False)
  SIM3.to_cbf("test_full_%s_003.cbf"%(params.context), intfile_scale=scale)
  diffs("CPU",SIM.raw_pixels, "GPU",SIM3.raw_pixels)

  print("\n# Use case 4 (%s): modularized api argchk=False, gpu_background=True"%(params.context))
  SIM4 = SWC.modularized_exafel_api_for_GPU(params=params, argchk=False, gpu_background=True)
  SIM4.to_cbf("test_full_%s_004.cbf"%(params.context), intfile_scale=scale)
  diffs("CPU",SIM.raw_pixels, "GPU",SIM4.raw_pixels)

  print("\n# Use case 5 (%s): modularized api argchk=True, gpu_background=True"%(params.context))
  SIM5 = SWC.modularized_exafel_api_for_GPU(params=params, argchk=True, gpu_background=True)
  SIM5.to_cbf("test_full_%s_005.cbf"%(params.context), intfile_scale=scale)
  diffs("CPU",SIM.raw_pixels, "GPU",SIM5.raw_pixels)

print("OK")
