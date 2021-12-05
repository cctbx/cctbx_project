from __future__ import absolute_import, division, print_function
from simtbx.nanoBragg import nanoBragg
import numpy as np

from simtbx.nanoBragg.tst_multisource_background import run_background_simulation

class gpu_run_background_simulation(run_background_simulation):

  def gpu_background(self, override_source=2):
    from simtbx.gpu import gpu_energy_channels
    gpu_channels_singleton = gpu_energy_channels(deviceId = 0)
    self.SIM.device_Id = 0
    # allocate GPU arrays
    from simtbx.gpu import exascale_api
    gpu_simulation = exascale_api(nanoBragg = self.SIM)
    gpu_simulation.allocate_cuda()

    from simtbx.gpu import gpu_detector as gpud
    gpu_detector = gpud(deviceId=self.SIM.device_Id, nanoBragg=self.SIM)
    gpu_detector.each_image_allocate()

    per_image_scale_factor = 0.0
    gpu_detector.scale_in_place(per_image_scale_factor) # apply scale directly on GPU
    gpu_simulation.add_background_cuda(gpu_detector)
    gpu_detector.write_raw_pixels(self.SIM)  # updates SIM.raw_pixels from GPU
    self.bg_multi = self.SIM.raw_pixels.as_numpy_array()

    gpu_detector.scale_in_place(per_image_scale_factor)
    gpu_simulation.add_background_cuda(detector = gpu_detector, override_source=override_source)
    gpu_detector.write_raw_pixels(self.SIM)
    self.bg_single = self.SIM.raw_pixels.as_numpy_array()

    gpu_detector.each_image_free()

def c_g_validate(cpu,gpu):
  mean_c = cpu.mean()
  mean_g = gpu.mean()
  print("cpu mean: %1.5g" % mean_c)
  print("gpu mean: %1.5g" % mean_g)
  if np.allclose(mean_c, mean_g): return True
  else:
    frac = mean_c / mean_g
    print("Means are off by a factor of %.6f" % frac)
    return False

def test_cpu_gpu_equivalence(n_chan=5,wave_interval=(0.998, 1.002),spectrum='tophat',single=2):
  cpu_run = run_background_simulation()
  gpu_run = gpu_run_background_simulation()
  beam = cpu_run.make_multichannel_beam_simulation(n_chan,wave_interval,spectrum)
  cpu_run.set_beam(beam)
  gpu_run.set_beam(beam)
  cpu_run.cpu_background(override_source=single)
  gpu_run.gpu_background(override_source=single)
  assert c_g_validate(cpu_run.bg_single, gpu_run.bg_single)
  if "plot" in sys.argv: plot_one_and_multi(cpu_run.bg_single, gpu_run.bg_single)
  assert c_g_validate(cpu_run.bg_multi, gpu_run.bg_multi)
  if "plot" in sys.argv: plot_one_and_multi(cpu_run.bg_multi, gpu_run.bg_multi)

if __name__=="__main__":
  import sys

  print("test with thin bandpass and tophat spectrum")
  test_cpu_gpu_equivalence()

  print("test with thin bandpass and tophat spectrum, more channels")
  test_cpu_gpu_equivalence(n_chan = 10)

  print("test with wider bandpass and tophat spectrum")
  test_cpu_gpu_equivalence(wave_interval=(0.98, 1.02))

  print("test with thin bandpass and gaussian spectrum")
  test_cpu_gpu_equivalence(n_chan=5, spectrum='gaussian', single=2)

  print("Scale-up the GPU channels--prints timing")
  from libtbx.development.timers import Profiler
  for n_chan in [20,40,80,160,320,640,1280,2560,5120]:
    P = Profiler("%d channels"%n_chan)
    gpu_run = gpu_run_background_simulation()
    beam = gpu_run.make_multichannel_beam_simulation(n_chan,wave_interval=(0.90,1.1))
    gpu_run.set_beam(beam)
    gpu_run.gpu_background()
    del P

  print("OK")
