from __future__ import absolute_import, division, print_function
from collections import namedtuple

from cctbx.xray.structure_factors.from_scatterers_direct \
  import from_scatterers_direct
from cctbx.xray.structure_factors.from_scatterers_fft \
  import from_scatterers_fft
from cctbx.xray.structure_factors.gradients_direct \
  import gradients_direct
from cctbx.xray.structure_factors.gradients_fft \
  import gradients_fft


Algorithm = namedtuple(
    "Algorithm",
    ["name", "desc", "from_scatterers", "gradients"]
)


algorithms = {
    "direct" : Algorithm(
        "direct",
        "Direct Summation",
        from_scatterers_direct,
        gradients_direct
    ),
    "fft" : Algorithm(
        "fft",
        "FFT Approximation",
        from_scatterers_fft,
        gradients_fft
    )
}

# discamb
try:
    from pydiscamb.cctbx_interface import from_scatterers_taam, gradients_taam
    algorithms["taam"] = Algorithm(
        "taam",
        "Transferrable Aspherical Atom Model",
        from_scatterers_taam,
        gradients_taam
    )
except ImportError:
    pass
