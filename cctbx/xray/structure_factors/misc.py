from __future__ import absolute_import, division, print_function
from cctbx import adptbx

def quality_factor_from_any(d_min=None,
                            grid_resolution_factor=None,
                            quality_factor=None,
                            u_base=None,
                            b_base=None):
  assert [quality_factor, u_base, b_base].count(None) >= 2
  if (u_base is not None):
    b_base = adptbx.u_as_b(u_base)
  if (b_base is not None):
    assert [d_min, grid_resolution_factor].count(None) == 0
    assert d_min > 0
    sigma = 1 / (2. * grid_resolution_factor)
    log_quality_factor = b_base * sigma * (sigma - 1) / (d_min * d_min)
    quality_factor = 10**log_quality_factor
  elif (quality_factor is None):
    quality_factor = 100
  return quality_factor

expensive_function_call_message = (
    "Programming problem: Calling this function is expensive."
  + " Please assign the result to an intermediate variable.")


# Check if pydiscamb is installed
try:
  # prefix to avoid including in structure_factors import
  import pydiscamb as __pydiscamb
  pydiscamb_is_installed = True
  del __pydiscamb
except ImportError:
  pydiscamb_is_installed = False
