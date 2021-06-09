"""
Verify the port from CUDA to KOKKOS
"""
from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from scitbx.matrix import sqr
from simtbx.nanoBragg import nanoBragg, shapetype

from simtbx.nanoBragg.tst_gauss_argchk import water, basic_crystal, basic_beam, basic_detector, amplitudes

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

 def several_wavelength_case_for_CPU(self):
  SIM = nanoBragg(self.DETECTOR, self.BEAM, panel_id=0)
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
  return SIM

 def modularized_exafel_api_for_KOKKOS(self, kokkos_run, argchk=False, cuda_background=True):
  from simtbx.kokkos import kokkos_energy_channels
  kokkos_channels_singleton = kokkos_energy_channels()

  SIM = nanoBragg(self.DETECTOR, self.BEAM, panel_id=0)
  SIM.device_Id = 0

  assert kokkos_channels_singleton.get_nchannels() == 0 # uninitialized
  assert kokkos_run.get_deviceID()==SIM.device_Id
  for x in range(len(self.flux)):
          kokkos_channels_singleton.structure_factors_to_KOKKOS_direct_cuda(
           x, self.sfall_channels[x].indices(), self.sfall_channels[x].data())
  assert kokkos_channels_singleton.get_nchannels() == len(self.flux)
  SIM.Ncells_abc = (20,20,20)
  SIM.Amatrix = sqr(self.CRYSTAL.get_A()).transpose()
  SIM.oversample = 2
  if argchk:
    print("\npolychromatic KOKKOS argchk")
    SIM.xtal_shape = shapetype.Gauss_argchk
  else:
    print("\npolychromatic KOKKOS no argchk")
    SIM.xtal_shape = shapetype.Gauss
  SIM.interpolate = 0


def kokkos_verification():
  from simtbx.kokkos import kokkos_instance, kokkos_energy_channels
  kokkos_run = kokkos_instance(deviceId = 0)

  kokkos_run.finalize_kokkos()
  print("RUN KOKKOS TRIAL")

if __name__=="__main__":
  # make the dxtbx objects
  BEAM = basic_beam()
  DETECTOR = basic_detector()
  CRYSTAL = basic_crystal()
  SF_model = amplitudes(CRYSTAL)
  # Famp = SF_model.Famp # simple uniform amplitudes
  SF_model.random_structure(CRYSTAL)
  SF_model.ersatz_correct_to_P1()

  #print("\n# Use case 2.  Three-wavelength polychromatic source")
  SWC = several_wavelength_case(BEAM, DETECTOR, CRYSTAL, SF_model)
  SIM = SWC.several_wavelength_case_for_CPU()
  #SIM.to_smv_format(fileout="test_full_e_002.img")
  #SIM.to_cbf("test_full_e_002.cbf")

  #print("\n# Use case: modularized api argchk=False, cuda_background=False")
  from simtbx.kokkos import kokkos_instance
  kokkos_run = kokkos_instance(deviceId = 0)
  SIM3 = SWC.modularized_exafel_api_for_KOKKOS(kokkos_run, argchk=False, cuda_background=False)
  kokkos_run.finalize_kokkos()
  #SIM3.to_cbf("test_full_e_003.cbf")
  #diffs("CPU",SIM.raw_pixels, "GPU",SIM3.raw_pixels) 

  print("OK")
