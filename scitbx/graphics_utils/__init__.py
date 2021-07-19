from __future__ import absolute_import, division, print_function
import boost_adaptbx.boost.python as bp
ext = bp.import_ext("scitbx_graphics_utils_ext")
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
    properties=properties.as_double(),
    selection=selection,
    color_all=color_all,
    gradient_type=gradients.index(gradient_type),
    min_value=min_value)


def map_to_rgb_colourmap(
    data_for_colors,
    colormap,
    selection,
    attenuation=None,
    powscale=1.0,
    map_directly=False,
    color_all=False):
  if attenuation is None:
    attenuation = flex.double(data_for_colors.size(), 1.0)
  return map_to_rgb_colourmap_(
    data_for_colors=data_for_colors.as_double(),
    colourmap=colormap,
    selection=selection,
    attenuation=attenuation,
    powscale=powscale,
    map_directly = map_directly,
    color_all=color_all)


def colour_by_phi_FOM(phi, fom = None):
  if not fom:
    fom = flex.double(phi.size(), 1.0)
  return color_by_phi_fom_(phi, fom)



def NoNansArray(arr, d=0.0):
  return NoNans(arr, d)

def IsNansArray(arr):
  return IsNans(arr)

def NoNansVecArray(vecs, defx=0.0, defy=0.0, defz=0.0):
  return NoNansvec3(vecs, defx, defy, defz)

def IsNansVecArray(vecs):
  return IsNansvec3(vecs)
