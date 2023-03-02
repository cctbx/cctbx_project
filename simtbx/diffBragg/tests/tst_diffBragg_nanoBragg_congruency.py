from __future__ import division
##from simtbx.kokkos import gpu_instance
#kokkos_run = gpu_instance(deviceId = 0)
from simtbx.diffBragg.utils import find_diffBragg_instances
import numpy as np


def main():
    from simtbx.nanoBragg.sim_data import SimData

    S = SimData(use_default_crystal=True)
    S.instantiate_diffBragg()
    S.D.nopolar = True
    S.D.oversample = 3

    S.D.add_diffBragg_spots()
    S._add_background()

    diff_img = S.D.raw_pixels.as_numpy_array()

    S.D.raw_pixels *= 0
    S.D.add_nanoBragg_spots()
    S._add_background()

    nano_img = S.D.raw_pixels.as_numpy_array()

    assert np.allclose(diff_img, nano_img, atol=1e-9)
    for name in find_diffBragg_instances(globals()): del globals()[name]


if __name__ == "__main__":
    import sys
    if "--kokkos" in sys.argv:
        import os
        os.environ["DIFFBRAGG_USE_KOKKOS"]="1"
    from simtbx.diffBragg.device import DeviceWrapper
    with DeviceWrapper(0) as _:
        main()
    print ("OK!")
