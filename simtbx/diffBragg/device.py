from __future__ import division
import libtbx.load_env  # implicit import
DIFFBRAGG_HAS_KOKKOS = libtbx.env.build_options.enable_kokkos

if DIFFBRAGG_HAS_KOKKOS:
    # TODO fix simtbx.kokkos import ...
    # from simtbx.kokkos import gpu_instance
    from simtbx.diffBragg import initialize_kokkos, finalize_kokkos

import os
USE_KOKKOS_GPU = DIFFBRAGG_HAS_KOKKOS and "DIFFBRAGG_USE_CUDA" in os.environ


class DeviceWrapper:

    def __init__(self, gpu_id):
        self.gpu_id = gpu_id
        self.instance = None

    def __enter__(self):
        if DIFFBRAGG_HAS_KOKKOS:
            if USE_KOKKOS_GPU:
                initialize_kokkos(self.gpu_id)
        return self

    def __exit__(self, tp, val, tbk):
        # TODO check for diffbragg instances
        if USE_KOKKOS_GPU:
            finalize_kokkos()
