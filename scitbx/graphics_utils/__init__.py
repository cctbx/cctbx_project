from __future__ import absolute_import, division, print_function
import boost.python
ext = boost.python.import_ext("scitbx_graphics_utils_ext")
from scitbx_graphics_utils_ext import *
from cctbx.array_family import flex

def color_by_property(
    properties,
    selection,
    color_all=False,
    gradient_type="rainbow",
    min_value=0.2):
  gradients = ["rainbow", "redblue", "heatmap"]
  assert (gradient_type in gradients)
  return color_by_property_(
    properties=properties,
    selection=selection,
    color_all=color_all,
    gradient_type=gradients.index(gradient_type),
    min_value=min_value)


def colour_by_phi_FOM(phi, fom = None):
  if not fom:
    fom = flex.double(phi.size(), 1.0)
  return color_by_phi_fom_(phi, fom)
