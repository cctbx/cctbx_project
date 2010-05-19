import boost.python
ext = boost.python.import_ext("cctbx_masks_ext")
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
    # XXX use scattering dictionary and set_selected
    # XXX use monomer library definitions for radii
    unknown = []
    if table is None:
      self.table = {}
    else:
      self.table = table.copy()
    self.atom_radii = flex.double()
    for i_seq, scatterer in enumerate(xray_structure.scatterers()):
      try:
        radius = self.table.get(scatterer.element_symbol())
        if radius is None:
          radius = van_der_waals_radii.vdw.table[scatterer.element_symbol()]
          self.table[scatterer.element_symbol()] = radius
      except KeyError:
        unknown.append(scatterer.element_symbol())
      else:
        self.atom_radii.append(radius)
    if(len(unknown) > 0):
      raise RuntimeError("Atoms with unknown van der Waals radius: ",unknown)

  def show(self, log=None):
    if log is None: log = sys.stdout
    for symbol in self.table:
      print >> log, "%5s" %symbol,
    print >> log
    for radius in self.table.values():
      print >> log, "%5.2f" %radius,
    print >> log


class _flood_fill(boost.python.injector, flood_fill):

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
    print >> log, "gridding: (%i,%i,%i)" %self.gridding_n_real()
    n_voids = self.n_voids()
    n_grid_points = reduce(lambda x,y:x*y, self.gridding_n_real())
    grid_points_per_void = self.grid_points_per_void()
    centres_of_mass = self.centres_of_mass_frac()
    eigensystems = self.eigensystems_frac()
    if self.n_voids() == 0: return
    print >> log, "Void #Grid points Vol/A^3 Vol/%  Centre of mass (frac)",
    print >> log, "  Eigenvectors (frac)"
    for i in range(n_voids):
      void_vol = (
        self.unit_cell().volume() * grid_points_per_void[i]) / n_grid_points
      formatted_site = ["%6.3f" % x for x in centres_of_mass[i]]
      print >> log, "%4i" %(i+1),
      print >> log, "%12i" %(grid_points_per_void[i]),
      print >> log, "%7.1f" %(void_vol),
      print >> log, "%5.1f" %(100*void_vol/self.unit_cell().volume()),
      print >> log, " (%s)" %(','.join(formatted_site)),
      for j in range(3):
        formatted_eigenvector = [
          "%6.3f" % x for x in eigensystems[i].vectors()[3*j:3*j+3]]
        if j == 0: sep = ""
        else: sep = " "*56
        print >> log, sep, "%i" %(j+1),
        print >> log, " (%s)" %(','.join(formatted_eigenvector))
