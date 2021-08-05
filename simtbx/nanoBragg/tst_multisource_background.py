from __future__ import absolute_import, division, print_function
from simtbx.nanoBragg import nanoBragg, nanoBragg_beam
from dials.array_family import flex
import numpy as np
"""Purpose of the test:  compare nanoBragg background two ways:
1) single channel
2) multiple channels
Overall photon fluence is the same in both simulations.
Results will be nearly identical if the multiple channel bandpass is small,
and if the spectrum is even (tophat), not irregular (random).
"""

water = flex.vec2_double([(0,2.57),(0.0365,2.58),(0.07,2.8),(0.12,5),(0.162,8),(0.18,7.32),(0.2,6.75),(0.216,6.75),(0.236,6.5),(0.28,4.5),(0.3,4.3),(0.345,4.36),(0.436,3.77),(0.5,3.17)])

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

class run_background_simulation:
  def __init__(self):
    self.SIM = nanoBragg()
    self.SIM.progress_meter = False
    self.SIM.Fbg_vs_stol = water
    self.SIM.amorphous_sample_thick_mm = 0.1
    self.SIM.amorphous_density_gcm3 = 1
    self.SIM.amorphous_molecular_weight_Da = 18
    self.total_flux = self.SIM.flux = 1e12
    self.verbose_beams = False

  def make_multichannel_beam_simulation(self,
      n_chan=5, wave_interval=(0.998, 1.002), spectrum='tophat'):
    assert spectrum in ['tophat','random','gaussian']
    beam = nanoBragg_beam.NBbeam()
    waves = np.linspace(wave_interval[0], wave_interval[1], n_chan)
    if spectrum=='tophat': fluences = np.ones(n_chan)
    elif spectrum=='random': fluences = np.random.random(n_chan)
    else: mu=(n_chan-1)/2.; sig=(n_chan-1)/6.; fluences = np.array([
      gaussian(i,mu,sig) for i in range(n_chan)])
    fluences /= fluences.sum() # sum of values is 1.
    fluences *= self.total_flux # sum of values is SIM.flux
    assert np.allclose(fluences.sum(), self.total_flux)
    beam.spectrum = list(zip(waves, fluences))
    return beam

  def set_beam(self,beam):
    self.SIM.verbose= 10 if self.verbose_beams else 0
    self.SIM.xray_beams = beam.xray_beams
    self.SIM.verbose=0
    if beam._undo_nanoBragg_norm_by_nbeams:
      assert np.allclose(self.SIM.flux, self.total_flux / len(beam.xray_beams))
    else:
      assert np.allclose(self.SIM.flux, self.total_flux)

  def cpu_background(self,override_source=2):
    self.SIM.raw_pixels *= 0
    self.SIM.add_background()
    self.bg_multi = self.SIM.raw_pixels.as_numpy_array()

    self.SIM.raw_pixels *= 0
    self.SIM.add_background(oversample=1, override_source=override_source)
    self.bg_single = self.SIM.raw_pixels.as_numpy_array()

def validate(multi,single): # range of sources or single source
  mean_single = single.mean()
  mean_multi = multi.mean()
  print("single source mean: %1.5g" % mean_single)
  print("multi source mean: %1.5g" % mean_multi)
  if np.allclose(mean_single, mean_multi): return True
  else:
    frac = mean_multi / mean_multi
    print("Means are off by a factor of %.6f" % frac)
    return False

def plot_one_and_multi(multi,single):
  from matplotlib import pyplot as plt
  fig,ax = plt.subplots(1,3)
  scale = multi.max()
  im = ax[0].imshow(multi, vmin=-scale, vmax=scale); ax[0].set_title("All sources")
  ax[1].imshow(single, vmin=-scale, vmax=scale); ax[1].set_title("Single source")
  ax[2].imshow(multi-single, vmin=-scale, vmax=scale); ax[2].set_title("Difference")
  fig.subplots_adjust(right=0.88)
  cbar = fig.add_axes([0.90,0.2,0.04,0.6]) # left, bottom, width, height
  fig.colorbar(im, cax=cbar)
  plt.show()

if __name__=="__main__":
  import sys
  run1 = run_background_simulation()
  nbeam_norm_check = [] # all test cases should give approx equal background as flux is constant

  print("test with thin bandpass and tophat spectrum")
  beam = run1.make_multichannel_beam_simulation(n_chan=5)
  run1.set_beam(beam)
  run1.cpu_background()
  nbeam_norm_check.append(run1.bg_multi.mean())
  assert validate(run1.bg_multi, run1.bg_single)
  if "plot" in sys.argv: plot_one_and_multi(run1.bg_multi, run1.bg_single)

  print("test with thin bandpass and tophat spectrum, more channels -- should fail")
  beam = run1.make_multichannel_beam_simulation(n_chan=10)
  run1.set_beam(beam)
  run1.cpu_background()
  nbeam_norm_check.append(run1.bg_multi.mean())
  assert not validate(run1.bg_multi, run1.bg_single) # fail since single-source isn't central
  if "plot" in sys.argv: plot_one_and_multi(run1.bg_multi, run1.bg_single)

  print("test with wider bandpass and tophat spectrum -- should fail")
  beam = run1.make_multichannel_beam_simulation(wave_interval=(0.98, 1.02))
  run1.set_beam(beam)
  run1.cpu_background()
  nbeam_norm_check.append(run1.bg_multi.mean())
  assert not validate(run1.bg_multi, run1.bg_single)
  if "plot" in sys.argv: plot_one_and_multi(run1.bg_multi, run1.bg_single)

  print("test with thin bandpass and gaussian spectrum -- it works")
  beam = run1.make_multichannel_beam_simulation(n_chan=15, spectrum='gaussian')
  run1.set_beam(beam)
  run1.cpu_background(override_source=7)
  nbeam_norm_check.append(run1.bg_multi.mean())
  assert validate(run1.bg_multi, run1.bg_single)

  print("compare mean background between tests: should be approx equal",nbeam_norm_check)
  assert np.allclose(nbeam_norm_check, nbeam_norm_check[0], atol=1.0)

  print("OK")

