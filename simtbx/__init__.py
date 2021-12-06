from __future__ import division

def get_exascale(interface, context):
  if context == "kokkos_gpu":
    from simtbx.kokkos import gpu_instance, gpu_energy_channels, gpu_detector, exascale_api
  elif context == "cuda":
    from simtbx.gpu import gpu_instance, gpu_energy_channels, gpu_detector, exascale_api
  else: raise NotImplementedError(context)

  return dict(gpu_instance = gpu_instance, gpu_energy_channels = gpu_energy_channels,
              gpu_detector = gpu_detector, exascale_api = exascale_api)[interface]

