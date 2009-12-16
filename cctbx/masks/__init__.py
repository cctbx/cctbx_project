import boost.python
ext = boost.python.import_ext("cctbx_masks_ext")
from cctbx_masks_ext import *

from cctbx.array_family import flex
from cctbx.eltbx import van_der_waals_radii

def vdw_radii_from_xray_structure(xray_structure):
  # XXX use scattering dictionary and set_selected
  # XXX use monomer library definitions for radii
  unknown = []
  atom_radii = flex.double()
  for i_seq, scatterer in enumerate(xray_structure.scatterers()):
    try:
      atom_radii.append(
        van_der_waals_radii.vdw.table[scatterer.element_symbol()])
    except:
      unknown.append(scatterer.element_symbol())
  if(len(unknown) > 0):
    raise RuntimeError("Atoms with unknown van der Waals radius: ",unknown)
  return atom_radii
