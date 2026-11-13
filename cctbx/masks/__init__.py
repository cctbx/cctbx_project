from __future__ import absolute_import, division, print_function
import boost_adaptbx.boost.python as bp
from functools import reduce
from six.moves import range
ext = bp.import_ext("cctbx_masks_ext")
from cctbx_masks_ext import *

import sys

from cctbx.array_family import flex
from cctbx.eltbx import van_der_waals_radii
from scitbx import linalg

# This function is still used in mmtbx/masks.py
def vdw_radii_from_xray_structure(xray_structure, table=None):
  return vdw_radii(xray_structure, table).atom_radii

class vdw_radii:

  def __init__(self, xray_structure, table=None):
    # XXX use monomer library definitions for radii
    unknown = []
    if table is None:
      self.table = {}
    else:
      self.table = table.copy()
    from cctbx.xray import ext as xray_ext
    scatterers = xray_structure.scatterers()
    # one lookup per scattering type, which decides the element symbol; a
    # registry of its own because the structure's would assign tables
    registry = xray_ext.scattering_type_registry()
    registry.process(scatterers)
    unique_indices = registry.unique_indices(scatterers)
    self.atom_radii = flex.double(unique_indices.size(), 0)
    unknown_types = set()
    for scattering_type, i in registry.type_index_pairs_as_dict().items():
      selection = (unique_indices == i)
      symbol = scatterers[selection.iselection()[0]].element_symbol()
      radius = self.table.get(symbol)
      if radius is None:
        if symbol not in van_der_waals_radii.vdw.table:
          unknown_types.add(scattering_type)
          continue
        radius = van_der_waals_radii.vdw.table[symbol]
        self.table[symbol] = radius
      self.atom_radii.set_selected(selection, radius)
    if(len(unknown_types) > 0):
      for scatterer in scatterers:
        if scatterer.scattering_type in unknown_types:
          unknown.append(scatterer.element_symbol())
      raise RuntimeError("Atoms with unknown van der Waals radius: ",unknown)

  def show(self, log=None):
    if log is None: log = sys.stdout
    symbols = list(self.table.keys())
    symbols.sort()
    print(" ".join("%5s" %symbol for symbol in symbols), file=log)
    print(" ".join("%5.2f" %self.table[symbol] for symbol in symbols), file=log)


@bp.inject_into(flood_fill)
class _():

  def eigensystems_frac(self):
    inertia_tensors = self.inertia_tensors_frac()
    return [linalg.eigensystem_real_symmetric(tensor)
            for tensor in inertia_tensors]

  def eigensystems_cart(self):
    inertia_tensors = self.inertia_tensors_cart()
    return [linalg.eigensystem_real_symmetric(tensor)
            for tensor in inertia_tensors]

  def show_summary(self, log=None):
    if log is None: log = sys.stdout
    print("gridding: (%i,%i,%i)" %self.gridding_n_real(), file=log)
    n_voids = self.n_voids()
    n_grid_points = reduce(lambda x,y:x*y, self.gridding_n_real())
    grid_points_per_void = self.grid_points_per_void()
    centres_of_mass = self.centres_of_mass_frac()
    eigensystems = self.eigensystems_frac()
    if self.n_voids() == 0: return
    print("Void #Grid points Vol/A^3 Vol/%  Centre of mass (frac)", end=' ', file=log)
    print("  Eigenvectors (frac)", file=log)
    for i in range(n_voids):
      void_vol = (
        self.unit_cell().volume() * grid_points_per_void[i]) / n_grid_points
      formatted_site = ["%6.3f" % x for x in centres_of_mass[i]]
      print("%4i" %(i+1), end=' ', file=log)
      print("%12i" %(grid_points_per_void[i]), end=' ', file=log)
      print("%7.1f" %(void_vol), end=' ', file=log)
      print("%5.1f" %(100*void_vol/self.unit_cell().volume()), end=' ', file=log)
      print(" (%s)" %(','.join(formatted_site)), end=' ', file=log)
      for j in range(3):
        formatted_eigenvector = [
          "%6.3f" % x for x in eigensystems[i].vectors()[3*j:3*j+3]]
        if j == 0: sep = ""
        else: sep = " "*56
        print(sep, "%i" %(j+1), end=' ', file=log)
        print(" (%s)" %(','.join(formatted_eigenvector)), file=log)
