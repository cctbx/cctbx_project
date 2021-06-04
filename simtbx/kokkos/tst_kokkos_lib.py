"""
Verify the port from CUDA to KOKKOS
"""
from __future__ import absolute_import, division, print_function



def kokkos_verification():
  from simtbx.kokkos import kokkos_instance
  kokkos_run = kokkos_instance(deviceId = 0)

  kokkos_run.finalize_kokkos()
  print("RUN KOKKOS TRIAL")


if __name__=="__main__":
  print("\n# Use case: dry-run kokkos")
  kokkos_verification()

print("OK")
